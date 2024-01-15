from PySide6.QtCore import *
from PySide6.QtSql import *

class RNASuiteDataBackend(QObject):
	error = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)

	def create_tables(self):
		sql = """
		CREATE TABLE outputs (
			id INTEGER PRIMARY KEY,
			type TEXT,
			name TEXT,
			
		)
		"""

	def connect_to_db(self, file=':memory:'):
		self.db = QSqlDatabase.addDatabase("QSQLITE")
		self.db.setDatabaseName(file)

		if not self.db.open():
			self.error.emit("Database connection error: {}".format(
				self.db.lastError().databaseText()
			))
			return




