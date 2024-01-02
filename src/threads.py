import os

from rchitect.utils import Rhome
from PySide6.QtCore import *

__all__ = ['RNASuiteCranMirrorThread', 'RNASuitePackageVersionThread']

class RNASuiteRbaseThread(QThread):
	result = Signal(object)
	error = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)

	def on_error_occurred(self, msg):
		self.error.emit(str(msg))

	def get_rscript(self):
		if os.name == 'nt':
			self.rscript = os.path.join(Rhome(), 'bin', 'Rscript.exe')
		else:
			self.rscript = os.path.join(Rhome(), 'bin', 'Rscript')

		if os.path.isfile(self.rscript):
			return True
		else:
			self.on_error_occurred("Could not find Rscript: {}".format(self.rscript))
			return False

	def process(self):
		pass

	def run(self):
		if self.get_rscript():
			self.process()

class RNASuiteCranMirrorThread(RNASuiteRbaseThread):
	def process(self):
		rcode = (
			"df=getCRANmirrors();"
			"write.table(df, sep='\t', quote=FALSE)"
		)
		proc = QProcess()
		proc.errorOccurred.connect(self.on_error_occurred)
		proc.start(self.rscript, ['-e', rcode])
		proc.waitForFinished()
		res = proc.readAllStandardOutput()
		lines = res.data().decode().split('\n')

		if lines:
			self.result.emit(lines)

class RNASuitePackageVersionThread(RNASuiteRbaseThread):
	def process(self):
		rcode = (
			"rv<-list(name='R', ver=packageVersion('base'));"
			"write.table(rv, sep='\t', quote=F, row.names=F, col.names=F);"
			"ps<-installed.packages();"
			"ps<-ps[,c('Package','Version')];"
			"write.table(ps, sep='\t', quote=F, row.names=F, col.names=F)"
		)
		proc = QProcess()
		proc.errorOccurred.connect(self.on_error_occurred)
		proc.start(self.rscript, ['-e', rcode])
		proc.waitForFinished()
		res = proc.readAllStandardOutput()
		lines = res.data().decode()
		versions = {}

		for line in lines.split('\n'):
			if line.strip():
				p, v = line.split()
				versions[p.lower()] = v

		self.result.emit(versions)
