import pandas

from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from models import *
from config import *

__all__ = ['RNASuitePandasTable', 'RNASuiteTableWidgets']

class RNASuitePandasTable(QTableView):
	error = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)

		self._model = RNASuitePandasModel(self)
		self.setModel(self._model)

		self.setAlternatingRowColors(True)
		self.setCornerButtonEnabled(False)
		self.setSortingEnabled(True)
		self.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
		self.verticalHeader().setDefaultSectionSize(10)

	"""
	def read_file(self, data_file):
		if data_file.endswith('.csv'):
			self.data_frame = pandas.read_csv(data_file, index_col=0)

		elif data_file.endswith(('.tsv', '.txt')):
			self.data_frame = pandas.read_csv(data_file, index_col=0, sep='\t')

		elif data_file.endswith(('.xls', '.xlsx')):
			self.data_frame = pandas.read_excel(data_file, index_col=0)

		else:
			self.error.emit("Can not support this file format")
			self.data_frame = pandas.DataFrame()

		if not self.data_frame.empty:
			self.data_model.load_data(self.data_frame)
	"""

	def load_dataframe(self, df):
		self._model.load_data(df)

	def load_tight(self, tight):
		df = pandas.DataFrame.from_dict(tight, orient='tight')
		self._model.load_data(df)

	def load_data(self, data):
		if isinstance(pandas.DataFrame, data):
			self.load_dataframe(data)
		else:
			self.load_tight(data)

	def show_data(self, table, data_id):
		self._model.show_data(table, data_id)


class RNASuiteTableWidgets(QObject):
	intab = Signal(QWidget, str)
	outab = Signal(QWidget, str)
	title = Signal(QWidget, str)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.metas = RNASuiteDataTables()
		self.parent = parent
		self.tables = {}

	def has_table(self, name):
		return name in self.tables

	def is_empty(self, name):
		if self.has_table(name) and not self.tables[name].empty:
			return False

		return True

	def get_table(self, name, title=None):
		if name in self.tables:
			if title is not None:
				self.title.emit(self.tables[name], title)

			return self.tables[name]

		self.tables[name] = RNASuitePandasTable(self.parent)
		table = self.metas[name]

		if title is None:
			title = table.title

		if table.type == 1:
			self.intab.emit(self.tables[name], title)
		else:
			self.outab.emit(self.tables[name], title)

		return self.tables[name]

	def set_title(self, name, title):
		self.title.emit(self.tables[name], title)

	def get_tight(self, name):
		if name in self.tables:
			return self.tables[name].get_tight()

		return {}

	def get_data(self, name):
		if name in self.tables:
			return self.tables[name].get_data()

