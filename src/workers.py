import traceback

import pandas

from rchitect import *
from PySide6.QtCore import *

from utils import *

__all__ = ['RNASuiteEdgerIdentifyWorker', 'RNASuiteDeseqIdentifyWorker',
	'RNASuiteDeseqExtractWorker', 'RNASuiteEdgerExtractWorker',
	'RNASuiteDegsDistPlotWorker',
	'RNASuiteDEGVolcanoPlotWorker', 'RNASuiteDEGVennPlotWorker',
	'RNASuiteDEGUpsetPlotWorker'
]

class RNASuiteBaseWorker(QThread):
	error = Signal(str)

	script = None
	function = None
	rettype = None

	def __init__(self, parent=None, params=None):
		super().__init__(parent)
		self.parent = parent
		self.params = params

	def submit(self, data):
		self.parent.pyconn.send(data)

	def prepare_data(self):
		for k, v in self.params.items():
			self.params[k] = convert_dataframe_to_dict(v)

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
		try:
			self.source_file()
			self.prepare_data()
			self.submit({
				'action': 'call',
				'func': self.function,
				'params': self.params
			})

		except:
			self.error.emit(traceback.format_exc())
		

class RNASuiteDeseqIdentifyWorker(RNASuiteBaseWorker):
	script = 'R/deseq.R'
	function = 'rnasuite_deseq_find_degs'

class RNASuiteDeseqExtractWorker(RNASuiteBaseWorker):
	function = 'rnasuite_deseq_extract_degs'

class RNASuiteEdgerIdentifyWorker(RNASuiteBaseWorker):
	script = 'R/edger.R'
	function = 'rnasuite_edger_find_degs'

class RNASuiteEdgerExtractWorker(RNASuiteBaseWorker):
	function = 'rnasuite_edger_extract_degs'

class RNASuiteDegsDistPlotWorker(RNASuiteBaseWorker):
	script = 'R/distplot.R'
	function = 'rnasuite_degs_dist_plot_run'

class RNASuiteDEGVolcanoPlotWorker(RNASuiteBaseWorker):
	script = 'R/volcanoplot.R'
	function = 'rnasuite_degs_volcano_plot'

	@property
	def data(self):
		return {
			'volcanoplot_tool': self.params['tool'],
			'volcanoplot_top': self.params['top'],
			'volcanoplot_line': self.params['line'],
			'volcanoplot_gname': self.params['gname'],
			'volcanoplot_gnsep': self.params['gnsep'],
			'volcanoplot_gncol': self.params['gncol'],
			'volcanoplot_limit': self.params['limit'],
			'volcanoplot_colors': [self.params['ncolor'], self.params['ucolor'], self.params['dcolor']]
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'plot',
			'func': 'degs_volcano_plot'
		})

class RNASuiteDEGVennPlotWorker(RNASuiteBaseWorker):
	script = 'R/vennplot.R'
	rettype = 'plot'
	function = 'rnasuite_degs_venn_plot'

	@property
	def data(self):
		return {
			'vennplot_tool': self.params['tool'],
			'vennplot_contrasts': self.params['contrasts'],
			'vennplot_percent': self.params['percent'],
			'vennplot_colors': self.params['colors'],
			'vennplot_opacity': self.params['opacity'],
			'vennplot_degtype': self.params['degtype']
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'plot',
			'func': 'degs_venn_plot'
		})

class RNASuiteDEGUpsetPlotWorker(RNASuiteBaseWorker):
	script = 'R/upsetplot.R'
	rettype = 'plot'
	function = 'rnasuite_degs_upset_plot'

	@property
	def data(self):
		return {
			'upsetplot_tool': self.params['tool'],
			'upsetplot_contrasts': self.params['contrasts'],
			'upsetplot_degtype': self.params['degtype']
		}

	def run(self):
		self.submit({
			'action': 'call',
			'rtype': 'plot',
			'func': 'degs_upset_plot'
		})