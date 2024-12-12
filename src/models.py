import pandas

from PySide6.QtGui import *
from PySide6.QtCore import *

from utils import *
from backend import *

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
	_headers = []
	_fields = []
	_table = None

	def __init__(self, parent=None):
		super().__init__(parent)

		self._displays = []
		self._selected = []

		self._total_count = 0
		self._read_count = 0
		self._read_size = 200

		self._cache_data = {}

	def rowCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return len(self._displays)

	def columnCount(self, parent=QModelIndex()):
		if parent.isValid():
			return 0

		return len(self._headers)

	def data(self, index, role=Qt.DisplayRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			return self.get_value(row, col)

	def headerData(self, section, orientation, role=Qt.DisplayRole):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			return self._headers[section]

	def canFetchMore(self, parent):
		if parent.isValid():
			return False

		if self._read_count < self._total_count:
			return True

		return False

	def fetchMore(self, parent):
		if parent.isValid():
			return

		ids = RDB.get_column(self.read_sql)
		fetch_count = len(ids)
		fetch_end = self._read_count+fetch_count-1
		self.beginInsertRows(QModelIndex(), self._read_count, fetch_end)
		self._displays.extend(ids)
		self._read_count += fetch_count
		self.endInsertRows()

	@property
	def count_sql(self):
		return SqlQuery(self._table)\
			.select('COUNT(1)')\
			.first()

	@property
	def read_sql(self):
		remain_count = self._total_count - self._read_count
		fetch_count = min(self._read_size, remain_count)
		return SqlQuery(self._table)\
			.select('id')\
			.limit(fetch_count)\
			.offset(self._read_count)

	@property
	def get_sql(self):
		if self._fields:
			return SqlQuery(self._table)\
				.select(*self._fields)\
				.where('id=?')\
				.first()
		else:
			return SqlQuery(self._table)\
				.select()\
				.where('id=?')\
				.first()

	def __update_cache(self, row):
		row_id = self._displays[row]
		self._cache_data ={row: RDB.get_row(self.get_sql, row_id)}

	def __update_count(self):
		self.row_count.emit(self._total_count)
		self.col_count.emit(len(self._headers))

	def get_value(self, row, col):
		if row not in self._cache_data:
			self.__update_cache(row)

		return self._cache_data[row][col]

	def update(self):
		self.beginResetModel()
		self._read_count = 0
		self._selected = []
		self._total_count = RDB.get_one(self.count_sql)
		self._displays = RDB.get_column(self.read_sql)
		self._read_count = len(self._displays)
		self._cache_data = {}
		self.endResetModel()
		self.__update_count()

	def reset(self):
		self.beginResetModel()
		self._cache_data = {}
		self._read_count = 0
		self._displays = []
		self._selected = []
		self._total_count = 0
		self.endResetModel()
		self.__update_count()

	def clear(self):
		sql = SqlQuery(self.table_name).delete()
		RDB.query(sql)
		self.reset()

class RNASuiteOutputTreeModel(RNASuiteSqliteTableModel):
	_table = 'output'
	_headers = ['Name', '']
	_fields = ['name', 'status', 'type', 'dataid']

	def data(self, index, role=Qt.ItemDataRole):
		if not index.isValid():
			return None

		row = index.row()
		col = index.column()

		if role == Qt.DisplayRole:
			if col == 0:
				return self.get_value(row, col)

		elif role == Qt.DecorationRole:
			if col == 0:
				if self.get_value(row, 2):
					return QIcon('icons/chart.svg')

				else:
					return QIcon('icons/table.svg')

			elif col == 1:
				state = self.get_value(row, 1)

				if state == 1:
					return QIcon('icons/dot.svg')

				elif state == 2:
					return QIcon('icons/refresh.svg')

		return None

	@Slot(int)
	def update_status(self, rowid):
		state = self.get_value(rowid, 1)

		sql = SqlQuery(self._table)\
			.update(status=0)\
			.where('status=1')\
		RDB.query(sql)

		sql = SqlQuery(self._table)\
			.update(status=1)\
			.where('id=?')
		RDB.query(sql, rowid)

	@Slot(str)
	def update_refresh(self, name):
		sql = SqlQuery(self._table)\
			.update(status=2)\
			.where('name=?')
		RDB.query(sql, name)

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