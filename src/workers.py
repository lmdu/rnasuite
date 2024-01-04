import traceback

import pandas

from rchitect import *

from PySide6.QtCore import *

__all__ = ['RNASuiteEdgerDEGWorker', 'RNASuiteDeseqDEGWorker',
	'RNASuiteShowDEGWorker',
]

class RNASuiteBaseDEGWorker(QObject):
	error = Signal(str)
	finish = Signal()
	script = None

	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.data_mapping = {}

	def prepare_data(self):
		for var, val in self.data_mapping.items():
			self.parent.pyconn.send({
				'action': 'data',
				'dataframe': isinstance(val, dict),
				'variable': var,
				'value': val
			})

	def load_script(self):
		if self.script is None:
			return

		with open(self.script) as fs:
			code = fs.read()

		self.parent.pyconn.send({
			'action': 'eval',
			'code': code
		})

	def run(self):
		pass

	def start(self):
		try:
			self.load_script()
			self.prepare_data()
			self.run()
		except:
			self.error.emit(traceback.format_exc())
		finally:
			self.finish.emit()

class RNASuiteDeseqDEGWorker(RNASuiteBaseDEGWorker):
	script = 'R/deseq.R'

	def __init__(self, parent, read_counts, sample_info, params, gene_names=None):
		super().__init__(parent)
		self.parent = parent

		self.data_mapping = {
			'read_counts': read_counts,
			'sample_info': sample_info,
			'deseq_design': params['design'],
			'deseq_fdr': params['fdr'],
			'deseq_logfc': params['lgfc'],
			'deseq_contrast': [params['compare'], params['treatment'], params['control']]
		}

	def run(self):
		self.parent.pyconn.send({
			'action': 'call',
			'rtype': 'degs',
			'func': 'deseq_analysis_pipeline'
		})

class RNASuiteEdgerDEGWorker(RNASuiteBaseDEGWorker):
	script = 'R/edger.R'

	def __init__(self, parent, read_counts, sample_info, params, gene_names=None):
		super().__init__(parent)
		self.parent = parent

		self.data_mapping = {
			'read_counts': read_counts,
			'sample_info': sample_info,
			'edger_design': params['design'],
			'edger_fdr': params['fdr'],
			'edger_logfc': params['lgfc'],
			'edger_compare': params['compare'],
			'edger_replicate': params['replicate'],
			'edger_method': params['method'],
			'edger_bcv': params['bcv'],
			'edger_contrast': [params['treatment'], params['control']]
		}

	def run(self):
		self.parent.pyconn.send({
			'action': 'call',
			'rtype': 'degs',
			'func': 'edger_analysis_pipeline'
		})

class RNASuiteShowDEGWorker(RNASuiteBaseDEGWorker):
	def __init__(self, parent, params):
		super().__init__(parent)
		self.parent = parent
		self.tool = params['tool']

		if self.tool == 'deseq':
			self.data_mapping = {
				'deseq_fdr': params['fdr'],
				'deseq_logfc': params['lgfc'],
				'deseq_contrast': [
					params['compare'],
					params['treatment'],
					params['control']
				]
			}

		elif self.tool == 'edger':
			self.data_mapping = {
				'edger_fdr': params['fdr'],
				'edger_logfc': params['lgfc'],
				'edger_contrast': [
					params['treatment'],
					params['control']
				]
			}

	def run(self):
		func = '{}_show_degs'.format(self.tool)
		self.parent.pyconn.send({
			'action': 'call',
			'rtype': 'degs',
			'func': func
		})