import traceback

from PySide6.QtCore import *

from config import *

__all__ = ['RNASuiteErrorPrompter']

class RNASuiteErrorPrompter(QObject):
	error_occurred = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.errors = RNASuiteErrors()

	def prompt(self, name, *args):
		error = self.errors[name]

		if args:
			err_str = error.message.format(*args)
		else:
			err_str = error.message

		if error.traceback:
			err_str = "{}:\n{}".format(err_str, traceback.format_exc())

		self.error_occurred.emit(err_str)
