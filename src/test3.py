import sqlite3

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

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



	@property
	def count_sql(self):
		return "SELECT COUNT(1) FROM {} LIMIT 1".format(
			self.table_name
		)

	@property
	def read_sql(self):
		remain_count = self.total_count - self.read_count
		fetch_count = min(self.read_size, remain_count)
		return "SELECT id FROM {} LIMIT {},{}".format(
			self.table_name,
			self.read_count,
			fetch_count
		)

	@property
	def get_sql(self):
		return "SELECT * FROM {} WHERE id=? LIMIT 1".format(
			self.table_name
		)

	def update_cache(self, row):
		row_id = self.displays[row]
		self.cache_data ={row: DB.get_row(self.get_sql, row_id)}

	def get_value(self, row, col):
		if row not in self.cache_data:
			self.update_cache(row)

		return self.cache_data[row][col]

	def get_total(self):
		return DB.get_one(self.count_sql)

	def update(self):
		self.beginResetModel()
		self.read_count = 0
		self.selected = []
		self.total_count = self.get_total()
		self.displays = DB.get_column(self.read_sql)
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
		DB.delete()
	


