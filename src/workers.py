import traceback

import pandas

from rchitect import *

from PySide6.QtCore import *

__all__ = ['RNASuiteEdgerDEGWorker', 'RNASuiteDeseqDEGWorker',
	'RNASuiteShowDEGWorker', 'RNASuiteDEGDistPlotWorker'
]

class RNASuiteBaseDEGWorker(QObject):
	error = Signal(str)
	started = Signal()
	finished = Signal()
	script = None

	def __init__(self, parent=None, params=None):
		super().__init__(parent)
		self.parent = parent
		self.params = params

	@property
	def data(self):
		pass

	def submit(self, data):
		self.parent.pyconn.send(data)

	def prepare_data(self):
		for var, val in self.data.items():
			self.submit({
				'action': 'data',
				'dataframe': isinstance(val, dict),
				'variable': var,
				'value': val
			})

	def source_file(self):
		if self.script is None:
			return

		#file = QFile(self.script)
		#file.open(QIODevice.ReadOnly)
		#byte = file.readAll()
		#code = byte.data().decode()
		with open(self.script) as fh:
			code = fh.read()

		self.submit({
			'action': 'eval',
			'code': code
		})

	def run(self):
		pass

	def start(self):
		try:
			self.started.emit()
			self.source_file()
			self.prepare_data()
			self.run()

		except:
			self.error.emit(traceback.format_exc())

		finally:
			self.finished.emit()

class RNASuiteDeseqDEGWorker(RNASuiteBaseDEGWorker):
	script = 'R/deseq.R'

	def __init__(self, parent, read_counts, sample_info, params, gene_names=None):
		super().__init__(parent, params)
		self.read_counts = read_counts
		self.sample_info = sample_info
		self.gene_names = gene_names

	@property
	def data(self):
		return {
			'read_counts': self.read_counts,
			'sample_info': self.sample_info,
			'deseq_design': self.params['design'],
			'deseq_fdr': self.params['fdr'],
			'deseq_logfc': self.params['lgfc'],
			'deseq_contrast': [self.params['compare'], self.params['treatment'], self.params['control']]
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'degs',
			'func': 'deseq_analysis_pipeline'
		})

class RNASuiteEdgerDEGWorker(RNASuiteBaseDEGWorker):
	script = 'R/edger.R'

	def __init__(self, parent, read_counts, sample_info, params, gene_names=None):
		super().__init__(parent, params)
		self.read_counts = read_counts
		self.sample_info = sample_info
		self.gene_names = gene_names

	@property
	def data(self):
		return {
			'read_counts': self.read_counts,
			'sample_info': self.sample_info,
			'edger_design': self.params['design'],
			'edger_fdr': self.params['fdr'],
			'edger_logfc': self.params['lgfc'],
			'edger_compare': self.params['compare'],
			'edger_replicate': self.params['replicate'],
			'edger_method': self.params['method'],
			'edger_bcv': self.params['bcv'],
			'edger_contrast': [self.params['treatment'], self.params['control']]
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'degs',
			'func': 'edger_analysis_pipeline'
		})

class RNASuiteShowDEGWorker(RNASuiteBaseDEGWorker):
	@property
	def data(self):
		if self.params['tool'] == 'deseq':
			return {
				'deseq_fdr': self.params['fdr'],
				'deseq_logfc': self.params['lgfc'],
				'deseq_contrast': [
					self.params['compare'],
					self.params['treatment'],
					self.params['control']
				]
			}

		elif self.params['tool'] == 'edger':
			return {
				'edger_fdr': self.params['fdr'],
				'edger_logfc': self.params['lgfc'],
				'edger_contrast': [
					self.params['treatment'],
					self.params['control']
				]
			}

	def run(self):
		func = '{}_show_degs'.format(self.tool)
		self.submit({
			'action': 'call',
			'rtype': 'degs',
			'func': func
		})

class RNASuiteDEGDistPlotWorker(RNASuiteBaseDEGWorker):
	script = 'R/distplot.R'

	@property
	def data(self):
		return {
			'distplot_tool': self.params['tool'],
			'distplot_contrasts': self.params['contrasts'],
			'distplot_type': self.params['plot'],
			'distplot_label': self.params['showval'],
			'distplot_rotate': self.params['rotate'],
			'distplot_colors': [self.params['ucolor'], self.params['dcolor']]
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'plot',
			'func': 'degs_dist_plot'
		})
