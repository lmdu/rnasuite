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

	def load_data(self, data):
		self.beginResetModel()
		#self.display_count = 10
		self._data = data
		self.endResetModel()

class RNASuiteSqliteTableModel(QAbstractTableModel):
	row_count = Signal(int)
	col_count = Signal(int)
	header_name = []
	table_name = None

	def __init__(self, parent=None):
		super().__init__(parent)

		self.displays = []
		self.selected = []
		
		self.total_count = 0
		self.read_count = 0
		self.read_size = 200

		self.cache_data = {}

	def rowCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return len(self.displays)

	def columnCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		count = len(self.header_name)
		self.row_count.emit(count)
		return count

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			return self.datasets[row][col]

	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			return self.header_title[section]

	def canFetchMore(self, parent):
		if parent.isValid():
			return False

		if self.read_count < self.total_count:
			return True

		return False

	def fetchMore(self, parent):
		if parent.isValid():
			return

		ids = DB.get_column(self.read_sql)
		fetch_count = len(ids)
		fetch_end = self.read_count+fetch_count-1
		self.beginInsertRows(QModelIndex(), self.read_count, fetch_end)
		self.displays.extend(ids)
		self.read_count += fetch_count
		self.endInsertRows()

	@property
	def count_sql(self):
		return SqlQuery(self.table_name)\
			.select('COUNT(1)')\
			.first()

	@property
	def read_sql(self):
		remain_count = self.total_count - self.read_count
		fetch_count = min(self.read_size, remain_count)
		return SqlQuery(self.table_name)\
			.select('id')\
			.limit(fetch_count)\
			.offset(self.read_count)

	@property
	def get_sql(self):
		return SqlQuery(self.table_name)\
			.select('id')\
			.first()

	def update_cache(self, row):
		row_id = self.displays[row]
		self.cache_data ={row: RDB.get_row(self.get_sql, row_id)}

	def get_value(self, row, col):
		if row not in self.cache_data:
			self.update_cache(row)

		return self.cache_data[row][col]

	def update(self):
		self.beginResetModel()
		self.read_count = 0
		self.selected = []
		self.total_count = RDB.get_one(self.count_sql)
		self.displays = RDB.get_column(self.read_sql)
		self.read_count = len(self.displays)
		self.cache_data = {}
		self.endResetModel()
		self.row_count.emit(self.total_count)

	def reset(self):
		self.beginResetModel()
		self.cache_data = {}
		self.read_count = 0
		self.displays = []
		self.selected = []
		self.total_count = 0
		self.endResetModel()

	def clear(self):
		sql = SqlQuery(self.table_name).delete()
		RDB.query(sql)
		self.reset()


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