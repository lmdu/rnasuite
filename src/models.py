import json
import pandas

from PySide6.QtGui import *
from PySide6.QtCore import *

from utils import *
from backend import *

__all__ = [
	'RNASuitePandasModel',
	'RNASuiteInputTreeModel',
	'RNASuiteOutputTreeModel',
	'RNASuitePluginModel'
]

class RNASuitePluginModel(QAbstractTableModel):
	_headers = ['Name', 'Description', 'Version', 'Action']

	def __init__(self, parent=None):
		super().__init__(parent)
		self.datasets = []

	def columnCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			if self.datasets:
				return len(self.datasets[0])

		return 0

	def rowCount(self, parent=QModelIndex()):
		if parent == QModelIndex():
			return len(self.datasets)

		return 0

	def data(self, index, role=Qt.ItemDataRole):
		if not index.isValid():
			return None

		if role == Qt.DisplayRole:
			return self.datasets[index.row()][index.column()]

	def headerData(self, section, orientation, role):
		if role == Qt.DisplayRole:
			if orientation == Qt.Horizontal:
				return self._headers[section]

	def add_plugin(self, *args):
		self.beginResetModel()
		self.datasets.append(args)
		self.endResetModel()

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

	def show_data(self, table, data_id):
		sql = SqlQuery(table)\
			.select('data')\
			.where('id=?')\
			.first()

		data = RDB.get_one(sql, data_id)
		data = json.loads(data)
		data_frame = pandas.DataFrame.from_dict(data, orient='tight')
		self.load_data(data_frame)

class RNASuiteSqliteTableModel(QAbstractTableModel):
	row_count = Signal(int)
	col_count = Signal(int)
	_headers = []
	_fields = []
	_table = None

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

		if self.read_count < self.total_count:
			return True

		return False

	def fetchMore(self, parent):
		if parent.isValid():
			return

		ids = RDB.get_column(self.read_sql)
		fetch_count = len(ids)
		fetch_end = self.read_count+fetch_count-1
		self.beginInsertRows(QModelIndex(), self.read_count, fetch_end)
		self.displays.extend(ids)
		self.read_count += fetch_count
		self.endInsertRows()

	@property
	def count_sql(self):
		return SqlQuery(self._table)\
			.select('COUNT(1)')\
			.first()

	@property
	def read_sql(self):
		remain_count = self.total_count - self.read_count
		fetch_count = min(self.read_size, remain_count)
		return SqlQuery(self._table)\
			.select('id')\
			.limit(fetch_count)\
			.offset(self.read_count)

	@property
	def get_sql(self):
		#if self._fields:
		return SqlQuery(self._table)\
			.select(*self._fields)\
			.where('id=?')\
			.first()
		#else:
		#	return SqlQuery(self._table)\
		#		.select()\
		#		.where('id=?')\
		#		.first()

	def update_cache(self, row):
		row_id = self.displays[row]
		self.cache_data ={row: RDB.get_row(self.get_sql, row_id)}

	def update_count(self):
		self.row_count.emit(self.total_count)
		self.col_count.emit(len(self._headers))

	def get_value(self, row, col):
		if row not in self.cache_data:
			self.update_cache(row)

		return self.cache_data[row][col]

	def get_data_id(self, index):
		return self.displays[index.row()]

	def get_table(self):
		return self._table

	def update(self):
		self.beginResetModel()
		self.read_count = 0
		self.selected = []
		self.total_count = RDB.get_one(self.count_sql)
		self.displays = RDB.get_column(self.read_sql)
		self.read_count = len(self.displays)
		self.cache_data = {}
		self.endResetModel()
		self.update_count()

	def reset(self):
		self.beginResetModel()
		self.cache_data = {}
		self.read_count = 0
		self.displays = []
		self.selected = []
		self.total_count = 0
		self.endResetModel()
		self.update_count()

	def clear(self):
		sql = SqlQuery(self.table_name).delete()
		RDB.query(sql)
		self.reset()

	def update_row_by_dataid(self, data_id):
		if data_id in self.displays:
			row_id = self.displays.index(data_id)
			self.update_cache(row_id)

class RNASuiteOutputTreeModel(RNASuiteSqliteTableModel):
	_table = 'outputchart'
	_headers = ['Name', '']
	_fields = ['name', 'state', 'type']

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
				if self.get_value(row, 2) == 'plot':
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

	def get_output_by_name(self, name):
		sql = SqlQuery(self._table)\
			.select('id', 'type', 'code')\
			.where('name=?')\
			.first()
		return RDB.get_dict(sql, name)

	def update_output(self, result, outid):
		sql = SqlQuery(self._table)\
			.update('state', 'code', 'data')\
			.where('id=?')

		if result['type'] == 'plot':
			state = 1
			code = result['data']
			data = ''
		else:
			state = 2
			code = 0
			data = json.dumps(result['data'])

		RDB.update_row(sql, state, code, data, outid)
		self.update_row_by_dataid(outid)

	def add_output(self, result):
		if result['type'] == 'plot':
			state = 1
			code = result['data']
			data = ''
		else:
			state = 2
			code = 0
			data = json.dumps(result['data'])

		sql = SqlQuery(self._table)\
			.insert(7)
		RDB.insert_row(sql, None, result['name'], state,
			result['type'], result['kind'], code, data)

		"""
		sql = SqlQuery(self._table)\
			.select('id')\
			.where('name=?')\
			.first()
		data_id = RDB.get_one(sql, name)

		self.beginInsertRows(QModelIndex(), self.total_count-1, self.total_count)
		self.displays.append(data_id)
		self.total_count += 1
		self.endInsertRows()
		self.update_count()
		"""
		self.update()

	@Slot(int)
	def update_row_click_state(self, index):
		self.beginResetModel()
		sql = SqlQuery(self._table)\
			.update('state')\
			.where('state=?')
		RDB.update_row(sql, 1, 0)

		sql = SqlQuery(self._table)\
			.update('state')\
			.where('id=?')
		data_id = self.get_data_id(index)
		RDB.update_row(sql, 1, data_id)
		self.endResetModel()

class RNASuiteInputTreeModel(RNASuiteSqliteTableModel):
	_table = 'inputfile'
	_headers = ['Name', 'Tag']
	_fields = ['name', 'tag']

	def add_input(self, df, name, tag):
		rows = len(df)
		cols = len(df.columns)
		data = json.dumps(df.to_dict(orient='tight'))

		sql = SqlQuery(self._table)\
			.select('id')\
			.where('tag=?')\
			.first()
		data_id = RDB.get_one(sql, tag)

		if data_id:
			sql = SqlQuery(self._table)\
				.update('name', 'rows', 'cols', 'data')\
				.where('id=?')
			RDB.update_row(sql, name, rows, cols, data, data_id)
		else:
			sql = SqlQuery(self._table)\
				.insert(6)
			RDB.insert_row(sql, None, name, tag, rows, cols, data)

		self.update()

	def get_data(self, index):
		row_id = index.row()
		data_id = self.get_data_id(index)

		sql = SqlQuery(self._table)\
			.select('data')\
			.where('id=?')\
			.first()
		data = RDB.get_one(sql, data_id)
		return data

class RNASuiteOutputTreesModel(QAbstractTableModel):
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