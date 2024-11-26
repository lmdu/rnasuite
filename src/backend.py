import sqlite3

from PySide6.QtCore import *

from utils import *

__all__ = ['RDB', 'SqlQuery']

class DataBackend(QObject):
	error = Signal(str)
	mutex = QMutex()
	conn = None

	def __init__(self, parent=None):
		super().__init__(parent)

	def __del__(self):
		if self.conn is not None:
			self.conn.close()

	def connect(self, file=':memory:'):
		if self.conn is not None:
			self.conn.close()

		self.conn = sqlite3.connect(file)

	@property
	def cursor(self):
		self.mutex.lock()
		cur = self.conn.cursor()
		self.mutex.unlock()
		return cur

	def query(self, sql, args=None):
		if args:
			return self.cursor.execute(sql, args)
		else:
			return self.cursor.execute(sql)

	def get_one(self, sql, *args):
		for row in self.query(sql, args):
			return row[0]

	def get_row(self, sql, *args):
		for row in self.query(sql, args):
			return row

	def get_rows(self, sql, *args):
		for row in self.query(sql, args):
			yield row

	def get_column(self, sql, *args):
		return [row[0] for row in self.query(sql, args)]

	def get_dict(self, sql, *args):
		res = self.query(sql, args)
		fields = [col[0] for col in res.description]

		for row in res:
			return AttrDict(zip(fields, row))

	def get_dicts(self, sql, *args):
		res = self.query(sql, args)
		fields = [col[0] for col in res.description]

		for row in res:
			yield AttrDict(zip(fields, row))

	def get_fields(self, table):
		sql = "SELECT * FROM {} LIMIT 1".format(table)
		res = self.query(sql)
		return [col[0] for col in res.description]
	
	def insert(self, table, data):
		field_num = len(self.get_fields(table))
		sql = "INSERT INTO {} VALUES ({})".format(
			table,
			','.join(['?']*field_num)
		)
		self.cursor.executemany(sql, data)

	def delete(self, table, rowid):
		sql = "DELETE FROM {} WHERE id=?".format(table)
		self.query(sql, (rowid,))

	def update(self, table, rowid, **kwargs):
		keys = []
		vals = []

		for k, v in kwargs.items():
			keys.append("{}=?".format(k))
			vals.append(v)

		vals.append(rowid)
		sql = "UPDATE {} SET {} WHERE id=?".format(table, ','.join(keys))
		self.query(sql, vals)

class SqlQuery:
	def __init__(self, table):
		self._table = table
		self._selects = []
		self._updates = []
		self._inserts = []
		self._and_wheres = []
		self._or_wheres = []
		self._order_bys = []
		self._order_desc = False
		self._limit = []
		self._offset = []
		self._action = None
		self._querys = []

	def __add(self, item):
		if not self._querys:
			self._querys.append(self._action)

		self._querys.extend([' ', item])

	def select(self, *args):
		self._selects.extend(args)
		self._action = 'SELECT'
		return self

	def insert(self, *args):
		self._inserts.extend(args)
		self._action = 'INSERT'
		return self

	def update(self, *args):
		self._updates.extend(args)
		self._action = 'UPDATE'
		return self

	def delete(self):
		self._action = 'DELETE'
		return self

	def where(self, *args, logic='and'):
		if logic == 'and':
			self._and_wheres.extend(args)
		elif logic == 'or':
			self._or_wheres.extend(args)

		return self

	def orderby(self, *args, asc=True):
		self._order_bys.extend(args)
		self._order_asc = asc
		return self

	def limit(self, num):
		self._limit = num
		return self

	def offset(self, num):
		self._offset = num
		return self

	def build(self):
		match self._action:
			case 'SELECT':
				if self._selects:
					self.__add(','.join(self._selects))
				else:
					self.__add('*')

				self.__add("FROM {}".format(self._table))
			
			case 'INSERT':
				self.__add("INTO {}".format(self._table))

				if self._inserts and isinstance(self._inserts[0], int):
					self.__add("VALUES ({})".format(','.join(['?']*self._inserts[0])))
				else:
					self.__add("({}) VALUES ({})".format(
						','.join(self._inserts),
						','.join(['?']*len(self._inserts))
						)
					)

			case 'UPDATE':
				self.__add("{} SET".format(self._table))
				self.__add(','.join("{}=?".format(u) for u in self._updates))

			case 'DELETE':
				self.__add("FROM {}".format(self._table))

		if self._and_wheres or self._or_wheres:
			self.__add("WHERE")

			if self._and_wheres:
				self.__add(' AND '.join(self._and_wheres))

			if self._or_wheres:
				if self._and_wheres:
					self.__add('OR')

				self.__add(' OR '.join(self._or_wheres))

		if self._order_bys:
			if self._order_asc:
				self.__add("ORDER BY {}".format(','.join(self._order_bys)))
			else:
				self.__add("ORDER BY {} DESC".format(','.join(self._order_bys)))

		if self._limit:
			self.__add("LIMIT {}".format(self._limit))

		if self._offset:
			self.__add("OFFSET {}".format(self._offset))

		return ''.join(self._querys)


RDB = DataBackend()

