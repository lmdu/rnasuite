import traceback

import pandas

from rchitect import *

from PySide6.QtCore import *

__all__ = ['RNASuiteDEGWorker', 'RNASuiteShowDEGWorker',
]

class RNASuiteBaseWorker(QObject):
	error = Signal(str)
	finish = Signal()

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

	def run(self):
		pass

	def start(self):
		try:
			self.prepare_data()
			self.run()
		except:
			self.error.emit(traceback.format_exc())
		finally:
			self.finish.emit()

class RNASuiteDEGWorker(RNASuiteBaseWorker):
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

class RNASuiteShowDEGWorker(RNASuiteBaseWorker):
	def __init__(self, parent, params):
		super().__init__(parent)
		self.parent = parent

		self.data_mapping = {
			'deseq_fdr': params['fdr'],
			'deseq_logfc': params['lgfc'],
			'deseq_contrast': [params['compare'], params['treatment'], params['control']]
		}

	def run(self):
		self.parent.pyconn.send({
			'action': 'call',
			'rtype': 'degs',
			'func': 'deseq_show_degs'
		})