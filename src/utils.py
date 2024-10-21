import socket
import pandas

__all__ = ['AttrDict', 'ClassDict', 'get_unused_port', 'format_number_display',
	'RNASuiteError', 'RNASuiteTable', 'RNASuitePackage', 'compare_version',
	'convert_dataframe_to_dict', 'convert_dict_to_dataframe'
]

class ClassDict:
	def __getitem__(self, key):
		return getattr(self, key)

	def __iter__(self):
		for key in dir(self):
			if not key.startswith('_'):
				yield key

	def get_order(self):
		return getattr(self, '_orders')

class AttrDict(dict):
	def __getattr__(self, attr):
		return self[attr]

	def __setattr__(self, attr, val):
		self[attr] = val

class RNASuiteTable(AttrDict):
	def __init__(self, ttype, title):
		self.type = ttype
		self.title = title

class RNASuiteError(AttrDict):
	def __init__(self, message, traceback=False):
		self.message = message
		self.traceback = traceback

class RNASuitePackage(AttrDict):
	def __init__(self, name, description, repository, required):
		self.name = name
		self.description = description
		self.repository = repository
		self.required = required
		self.version = ''
		self.status = ''

def get_unused_port():
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	sock.bind(('', 0))
	_, port = sock.getsockname()
	sock.close()

	return port

def format_number_display(num):
	if isinstance(num, float):
		return '{:g}'.format(num)
	else:
		return str(num)

def compare_version(ver1, ver2):
	#ver1 is required version
	#ver2 is current version
	if not ver2:
		return 0

	if ver1 == ver2:
		return 1

	v1_num = [int(v) for v in ver1.split('.')]
	v2_num = [int(v) for v in ver2.split('.')]

	for i in range(max(len(v1_num), len(v2_num))):
		v1 = v1_num[i] if i < len(v1_num) else 0
		v2 = v2_num[i] if i < len(v2_num) else 0

		if v2 > v1:
			return 1
		elif v2 < v1:
			return -1

def convert_dataframe_to_dict(data):
	if isinstance(data, pandas.DataFrame):
		return data.to_dict(orient='tight')
	else:
		return data

def convert_dict_to_dataframe(data):
	if isinstance(data, dict):
		return pandas.DataFrame.from_dict(data, orient='tight')
	else:
		return data

if __name__ == '__main__':
	print(compare_version('2.0.2', '1.3.1'))
