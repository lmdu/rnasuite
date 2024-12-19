import os
import time
import functools
import traceback
import multiprocessing

import pandas
import rchitect
from rchitect.interface import set_hook, package_event, \
	peek_event, process_events, polled_events

from PySide6.QtCore import *

from utils import *
from config import *

__all__ = [
	'RNASuiteREnvironment',
	'RNASuiteRMessageProcessor'
]

class RNASuiteREnvironment(multiprocessing.Process):
	def __init__(self, rconn, port, token):
		super().__init__()
		self.name = "R session"
		self.rconn = rconn
		self.port = port
		self.token = token

	def send(self, action, content=None, data=None):
		self.rconn.send({
			'action': action,
			'content': content,
			'data': data
		})

	def source_file(self, rfile):
		root_dir = os.path.dirname(os.path.abspath(__file__))
		rfile = os.path.join(root_dir, 'R', rfile)
		rchitect.rcall('source', rfile)

	def set_hooks(self):
		hooks = {
			'base': 'rnasuite.R'
		}

		for package, script in hooks.items():
			set_hook(package_event(package, 'onLoad'), lambda *args: self.source_file(script))

	@staticmethod
	def console_callback(self, buf, otype):
		#if otype == 0:
		buf = buf.strip()
		if buf:
			self.send('message', buf)

	def set_callbacks(self):
		rchitect.def_callback(name='write_console_ex')(
			functools.partial(self.console_callback, self)
		)

	def get_package_version(self, package):
		try:
			return rchitect.rcopy(rchitect.rcall('packageVersion', package))
		except:
			return None

	def load_package(self, package):
		if self.get_package_version(package) is None:
			self.send('warning', "R package {} should be installed first".format(package))
			return False
		else:
			try:
				rchitect.rcall('library', package)
				return True
			except Exception:
				self.send('error', "Failed to load package {}:\n{}".format(package, traceback.format_exc()))
				return False

	def load_httpgd(self):
		if not self.load_package('httpgd'):
			return

		#start httpgd server
		try:
			rchitect.rcall('hgd',
				silent = False,
				port = self.port,
				token = self.token
			)
			self.send('socket')

		except Exception:
			self.send('error', "Failed to start httpgd:\n{}".format(traceback.format_exc()))

	def load_reticulate(self):
		if not self.load_package('reticulate'):
			return

		#self.source_file('pipe.R')

	def start_r_env(self):
		rchitect.init()
		#self.set_hooks()
		self.set_callbacks()

		self.load_httpgd()
		self.load_reticulate()
		self.source_file('rnasuite.R')
		#self.load_package('DESeq2')

	def run(self):
		self.start_r_env()

		while True:
			try:
				if self.rconn.poll():
					data = self.rconn.recv()

				else:
					if peek_event():
						process_events()

					else:
						polled_events()

					time.sleep(0.01)
					continue

			except (EOFError, BrokenPipeError):
				break

			try:
				self.send('running', data=True)

				if 'params' in data:
					for k, v in data['params'].items():
						data['params'][k] = convert_dict_to_dataframe(v)

				else:
					data['params'] = {}

				match data['action']:
					case 'data':
						if data['dataframe']:
							data['value'] = pandas.DataFrame.from_dict(data['value'], orient='tight')
							rchitect.rcall('send_df_to_r', data['variable'], data['value'])

						else:
							rchitect.rcall('send_val_to_r', data['variable'], data['value'])

					case 'call':
						ret = rchitect.rcopy(rchitect.rcall(data['func'], **data['params']))

						for i in range(len(ret)):
							if isinstance(ret[i]['data'], pandas.DataFrame):
								ret[i]['data'] = convert_dataframe_to_dict(ret[i]['data'])

						self.send('result', data=ret)

					case 'eval':
						rchitect.reval(data['code'])

					case 'plot':
						ret = rchitect.rcopy(rchitect.rcall('unigd::ugd_render', **data['params']))
						self.send('plot', ret)

					case 'remove':
						rchitect.rcall('unigd::ugd_remove', **data['params'])

			except:
				self.send('error', traceback.format_exc())

			finally:
				self.send('running', data=False)

class RNASuiteRMessageProcessor(QThread):
	socket = Signal()
	plot = Signal(str)
	error = Signal(str)
	warning = Signal(str)
	message = Signal(str)
	results = Signal(object)
	running = Signal(bool)

	def __init__(self, parent):
		super().__init__(parent)
		self.parent = parent

	def data_transfrom(self, data):
		
		if isinstance(data, dict):
			#if 'index_names' in data and 'column_names' in data:
			return pandas.DataFrame.from_dict(data, orient='tight')

		return data

	def run(self):
		while True:
			try:
				data = self.parent.pyconn.recv()

			except EOFError:
				break

			match data['action']:
				case 'warning':
					self.warning.emit(data['content'])

				case 'error':
					self.error.emit(data['content'])
				
				case 'message':
					self.message.emit(data['content'])

				case 'running':
					self.running.emit(data['data'])

				case 'plot':
					self.plot.emit(data['content'])

				case 'result':
					self.results.emit(data['data'])

				case 'socket':
					self.socket.emit()

