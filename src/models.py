import pandas

from PySide6.QtGui import *
from PySide6.QtCore import *

from utils import *

__all__ = ['RNASuitePandasModel', 'RNASuiteOutputTreeModel']

class RNASuitePandasModel(QAbstractTableModel):
	def __init__(self, parent=None):
		super().__init__(parent)
		self._data = pandas.DataFrame()

	def rowCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			return len(self._data)

		return 0

	def columnCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			return len(self._data.columns)

		return 0

	def data(self, index, role=Qt.ItemDataRole):
		if not index.isValid():
			return None

		if role == Qt.DisplayRole:
			v = self._data.iloc[index.row(), index.column()]
			return format_number_display(v)

		return None


	def headerData(self, section, orientation, role):
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				return self._data.columns[section]

			if orientation == Qt.Vertical:
				return self._data.index[section]

		return None

	'''
	def canFetchMore(self, parent=QModelIndex()):
		return self.display_count < len(self._data)

	def fetchMore(self, parent=QModelIndex()):
		add_count = 100

		if self.display_count + add_count > len(self._data):
			add_count = len(self._data) - self.display_count

		self.beginInsertRows(QModelIndex(), self.display_count, self.display_count+add_count-1)
		self.display_count += add_count
		self.endInsertRows()
	'''

	"""
	def sort(self, column, order):
		if self._data.empty:
			return

		if order ==  Qt.DescendingOrder:
			ascending = False
		else:
			ascending = True

		self.beginResetModel()
		self._data.sort_values(self._data.columns[column],
			ascending = ascending,
			inplace = True
		)
		self.endResetModel()
	"""

	def load_data(self, data):
		self.beginResetModel()
		#self.display_count = 10
		self._data = data
		self.endResetModel()

class RNASuiteOutputTreeModel(QAbstractTableModel):
	_headers = ['Name', '', 'plot', 'id', 'type', 'pyid']

	def __init__(self, parent=None):
		super().__init__(parent)
		self._data = IndexedDict()

	def add_row(name, status, _type):
		self._data[name] = AttrDict(
			name = name,
			type = _type,
			status = status,
		)

	def update_status(name, status):
		self._data[name].status = status

	def rowCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			return len(self._data)

		return 0

	def columnCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			return len(self._headers)

		return 0

	def data(self, index, role=Qt.ItemDataRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			if col == 0:
				v = self._data[row]
				return v.name

		elif role == Qt.DecorationRole:
			if col == 0:
				v = self._data[row]

				if v.type == 'plot':
					return QIcon('icons/chart.svg')

				else:
					return QIcon('icons/table.svg')

			elif col == 1:
				v = self._data[row]

				if v.type == 'table' and self.status:
					return QIcon('icons/refresh.svg')

				if v.type == 'plot' and self.status:
					return QIcon('icons/dot.svg')

		return None


	def headerData(self, section, orientation, role):
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				return self._headers[section]

			if orientation == Qt.Vertical:
				return self._data.index[section]

		return None