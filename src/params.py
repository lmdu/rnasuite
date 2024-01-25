from utils import *

__all__ = [
	'RNASuiteDeseqParameters',
	'RNASuiteExtractDegsParameters',
	'RNASuiteEdgerParameters',
	'RNASuiteDegsDistPlotParameters',
	'RNASuiteDegsDistPlotControlParameters',
	'RNASuiteDegsVennPlotParameters',
	'RNASuiteDegsVennPlotControlParameters',
	'RNASuiteDegsUpsetPlotParameters',
	'RNASuiteDegsUpsetPlotControlParameters',
	'RNASuiteDegsVolcanoPlotParameters',
	'RNASuiteDegsVolcanoPlotControlParameters',
	'RNASuiteDeseqMaPlotControlParameters',
]

class RNASuiteParameter(dict):
	def __init__(self, expose=True, default=None, index=False, decimal=1, options=[], **kwargs):
		self.update(kwargs)
		self['expose'] = expose
		self['default'] = default
		self['index'] = index
		self['options'] = options
		self['decimal'] = decimal

	def __getattr__(self, key):
		return self[key]

class RNASuiteParameters:
	def __init__(self, *args):
		self.params = args
		self.mapping = {p.key: i for i, p in enumerate(self.params)}

	def __getitem__(self, key):
		return self.params[self.mapping[key]]

	def __contains__(self, key):
		return key in self.mapping

	def __iter__(self):
		for p in self.params:
			yield p

RNASuiteDeseqParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'fdr',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		decimal = 5,
		display = 'FDR:',
		default = 0.05
	),
	RNASuiteParameter(
		key = 'logfc',
		type = 'int',
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	RNASuiteParameter(
		key = 'compare',
		type = 'list',
		display = 'Comparison between:'
	),
	RNASuiteParameter(
		key = 'control',
		type = 'list',
		display = 'Control group:'
	),
	RNASuiteParameter(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
	),
	RNASuiteParameter(
		key = 'eliminate',
		type = 'select',
		display = 'Consider effect of factors:',
		expose = False
	),
	RNASuiteParameter(
		key = 'design',
		type = 'str',
		display = "Model design formula:"
	),
	RNASuiteParameter(
		key = 'custom',
		type = 'bool',
		display = "Custom design formula:",
		expose = False
	)
)

RNASuiteEdgerParameters = RNASuiteParameters(
	AttrDict(
		key = 'fdr',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		display = 'FDR:',
		default = 0.05
	),
	AttrDict(
		key = 'lgfc',
		type = 'int',
		range = (0, 100),
		step = 1,
		display = 'log2FoldChange:',
		default = 1
	),
	AttrDict(
		key = 'compare',
		type = 'list',
		options = [],
		display = 'Comparison between:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'control',
		type = 'list',
		options = [],
		display = 'Control group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'eliminate',
		type = 'select',
		options = [],
		display = 'Considering effect of factors:',
		default = None,
		index = False
	),
	AttrDict(
		key = 'design',
		type = 'str',
		display = "Model design formula:",
		default = None
	),
	AttrDict(
		key = 'custom',
		type = 'bool',
		display = "Custom design formula:",
		default = False
	),
	AttrDict(
		key = 'replicate',
		type = 'list',
		options = ['Biological Replicates', 'No Biological Replicates'],
		display = "Use edgeR with:",
		default = None,
		index = True
	),
	AttrDict(
		key = 'bcv',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		display = "BCV value:",
		default = 0.4
	),
	AttrDict(
		key = 'method',
		type = 'list',
		options = [
			"Quasi-likelihood F-tests",
			"Likelihood ratio tests",
			"Classic exact test"
		],
		display = "Testing method:",
		default = None,
		index = True
	)
)

RNASuiteExtractDegsParameters = RNASuiteParameters(
	# AttrDict(
	# 	key = 'fdr',
	# 	type = 'float',
	# 	range = (0, 1),
	# 	step = 0.01,
	# 	display = 'FDR:',
	# 	default = 0.05
	# ),
	# AttrDict(
	# 	key = 'lgfc',
	# 	type = 'int',
	# 	range = (0, 100),
	# 	step = 1,
	# 	display = 'log2FoldChange:',
	# 	default = 1
	# ),
	RNASuiteParameter(
		key = 'control',
		type = 'list',
		options = [],
		display = 'Control group:'
	),
	RNASuiteParameter(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:'
	)
)

RNASuiteDegsDistPlotParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'contrasts',
		type = 'contrast',
		default = None,
		display = 'Contrasts:'
	)
)

RNASuiteDegsVolcanoPlotParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'control',
		type = 'list',
		display = 'Control group:'
	),
	RNASuiteParameter(
		key = 'treatment',
		type = 'list',
		options = [],
		display = 'Treatment group:',
	),
	RNASuiteParameter(
		key = 'gname',
		type = 'list',
		options = ('gene ID column', 'gene annotation'),
		display = 'Show gene names from:',
		default = None,
		index = True
	),
	RNASuiteParameter(
		key = 'gnsep',
		type = 'str',
		display = "Gene ID separator:",
		default = '|'
	),
	RNASuiteParameter(
		key = 'gncol',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 1,
		display = "Show gene names in column:"
	)
)

RNASuiteDegsVennPlotParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'contrasts',
		type = 'contrast',
		display = 'Contrasts:'
	),
	RNASuiteParameter(
		key = 'degtype',
		type = 'list',
		display = "DEG Type:",
		options = ('All DEGs', 'Up-regulated', 'Down-regulated'),
		index = True
	)
)

RNASuiteDegsUpsetPlotParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'contrasts',
		type = 'contrast',
		display = 'Contrasts:'
	),
	RNASuiteParameter(
		key = 'degtype',
		type = 'list',
		display = "DEG Type:",
		options = ('All DEGs', 'Up-regulated', 'Down-regulated'),
		index = True
	)
)

RNASuiteDeseqMaPlotControlParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'main',
		type = 'str',
		default = '',
		display = "Title:"
	),
	RNASuiteParameter(
		key = 'xlab',
		type = 'str',
		default = 'mean of normalized counts',
		display = "X label:"
	),
	RNASuiteParameter(
		key = 'ylim',
		type = 'limit',
		range = (-1000, 1000),
		step = 1,
		default = (0, 0),
		display = "Y limit:"
	),
	RNASuiteParameter(
		key = 'colNonSig',
		type = 'color',
		default = '#999999',
		display = "Non-significant point color:"
	),
	RNASuiteParameter(
		key = 'colSig',
		type = 'color',
		default = '#0000FF',
		display = "Significant point color:"
	),
	RNASuiteParameter(
		key = 'colLine',
		type = 'color',
		default = '#666666',
		display = "Horizontal line color:"
	),
)

RNASuiteDegsDistPlotControlParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'plot_type',
		type = 'list',
		options = ['Stacked bar plot', 'Dodged bar plot'],
		index = True,
		display = "Plot type:"
	),
	RNASuiteParameter(
		key = 'show_label',
		type = 'bool',
		default = False,
		display = "Show value labels:"
	),
	RNASuiteParameter(
		key = 'bar_colors',
		type = 'colors',
		default = ('#E74C3C', '#3498DB'),
		display = "Up and down fill colors:"
	),
	RNASuiteParameter(
		key = 'x_rotate',
		type = 'int',
		range = [0, 90],
		step = 5,
		display = "X labels rotate angle:"
	),
	RNASuiteParameter(
		key = 'plot_title',
		type = 'str',
		display = "Title:"
	),
	RNASuiteParameter(
		key = 'x_label',
		type = 'str',
		display = "X label:"
	),
	RNASuiteParameter(
		key = 'y_label',
		type = 'str',
		display = "Y label:",
		default = "Counts"
	),
	RNASuiteParameter(
		key = 'legend_title',
		type = 'str',
		display = 'Legend title:',
		default = 'DEGs'
	),
	RNASuiteParameter(
		key = 'theme_name',
		type = 'list',
		options = ['bw', 'classic', 'linedraw',
			'minimal', 'void', 'light', 'grey',
			'gray', 'dark'
		],
		display = 'Plot theme:'
	),
	RNASuiteParameter(
		key = 'base_size',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 11,
		display = 'Base size:'
	)
)

RNASuiteDegsVennPlotControlParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'show_percentage',
		type = 'bool',
		display = "Show percent values:",
		default = True
	),
	RNASuiteParameter(
		key = 'digits',
		type = 'int',
		range = (0, 5),
		step = 1,
		default = 1,
		display = 'Decimal digits:'
	),
	RNASuiteParameter(
		key = 'fill_color',
		type = 'colors',
		default = ['#E74C3C', '#3498DB', '#2ECC71', '#F1C40F'],
		display = 'Fill colors:'
	),
	RNASuiteParameter(
		key = 'fill_alpha',
		type = 'float',
		range = (0, 1),
		step = 0.1,
		default = 0.5,
		display = "Fill alpha:"
	),
	RNASuiteParameter(
		key = 'stroke_color',
		type = 'color',
		display = "Storke color:",
		default = '#000000'
	),
	RNASuiteParameter(
		key = 'stroke_alpha',
		type = 'float',
		range = (0, 1),
		step = 0.1,
		default = 1,
		display = "Storke alpha:"
	),
	RNASuiteParameter(
		key = 'stroke_size',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 1,
		display = "Storke size:"
	),
	RNASuiteParameter(
		key = 'stroke_linetype',
		type = 'list',
		options = ['solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash'],
		display = "Storke line type:"
	),
	RNASuiteParameter(
		key = 'set_name_color',
		type = 'color',
		default = '#000000',
		display = "Set name color:"
	),
	RNASuiteParameter(
		key = 'set_name_size',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 6,
		display = "Set name size:"
	),
	RNASuiteParameter(
		key = 'text_color',
		type = 'color',
		default = '#000000',
		display = "Label color:"
	),
	RNASuiteParameter(
		key = 'text_size',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 4,
		display = "Label size:"
	)
)

RNASuiteDegsUpsetPlotControlParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'point_size',
		type = 'float',
		range = (0, 30),
		step = 0.1,
		default = 2.5,
		display = 'Matrix point size:'
	),
	RNASuiteParameter(
		key = 'point_size',
		type = 'float',
		range = (0, 30),
		step = 0.1,
		default = 2.5,
		display = 'Matrix point size:'
	),
	RNASuiteParameter(
		key = 'line_width',
		type = 'float',
		range = (0, 10),
		step = 0.1,
		default = 0.8,
		display = "Matrix line width:"
	),
	RNASuiteParameter(
		key = 'text_scale',
		type = 'float',
		range = (0, 10),
		step = 0.1,
		default = 1.8,
		display = "Text size:"
	),
	RNASuiteParameter(
		key = 'order_by',
		type = 'list',
		options = ('freq', 'degree'),
		display = "Order by:"
	),
	RNASuiteParameter(
		key = 'order_desc',
		type = 'bool',
		default = True,
		display = "Descreasing:"
	),
	RNASuiteParameter(
		key = 'intersect_num',
		type = 'int',
		range = (0, 1000),
		step = 1,
		default = 100,
		display = "Number of intersections:"
	),
	RNASuiteParameter(
		key = 'main_color',
		type = 'color',
		default = '#3b3b3b',
		display = 'Main bar color:'
	),
	RNASuiteParameter(
		key = 'main_ylabel',
		type = 'str',
		default = "Intersection Size",
		display = "Main y label:"
	),
	RNASuiteParameter(
		key = 'matrix_color',
		type = 'color',
		default = '#3b3b3b',
		display = 'Matrix dot color:'
	),
	RNASuiteParameter(
		key = 'matrix_alpha',
		type = 'float',
		range = (0, 1),
		step = 0.01,
		decimal = 2,
		default = 0.5,
		display = "Matrix dot alpha:"
	),
	RNASuiteParameter(
		key = 'set_color',
		type = 'color',
		default = '#3b3b3b',
		display = 'Set bar color:'
	),
	RNASuiteParameter(
		key = 'set_xlabel',
		type = 'str',
		default = "Set Size",
		display = "Set x label:"
	),
	RNASuiteParameter(
		key = 'show_number',
		type = 'list',
		options = ('yes', 'no'),
		display = "Show numbers:"
	),
	RNASuiteParameter(
		key = 'number_angles',
		type = 'int',
		range = (0, 360),
		step = 10,
		default = 0,
		display = "Number angles:"
	),
	RNASuiteParameter(
		key = 'show_empty',
		type = 'list',
		options = ('off', 'on'),
		display = "Show empty sets:"
	)
)

RNASuiteDegsVolcanoPlotControlParameters = RNASuiteParameters(
	RNASuiteParameter(
		key = 'top',
		type = 'int',
		range = (0, 1000),
		step = 1,
		default = 10,
		display = "Show top significant:"
	),
	RNASuiteParameter(
		key = 'point_color',
		type = 'colors',
		default = ('#ff006e', '#3a86ff', '#e9ecef'),
		display = "Point colors for Up/Down/NS:"
	),
	RNASuiteParameter(
		key = 'point_size',
		type = 'float',
		range = (0, 10),
		decimal = 1,
		step = 0.1,
		default = 0.5,
		display = "Point size:"
	),
	RNASuiteParameter(
		key = 'show_vline',
		type = 'bool',
		default = True,
		display = "Show vertical line:"
	),
	RNASuiteParameter(
		key = 'vline_type',
		type = 'list',
		options = ['dashed', 'dotted', 'dotdash', 'longdash', 'twodash', 'solid'],
		default = 'dashed',
		display = "Vertical line type:"
	),
	RNASuiteParameter(
		key = 'vline_color',
		type = 'color',
		default = '#000000',
		display = "Vertical line color:"
	),
	RNASuiteParameter(
		key = 'vline_width',
		type = 'float',
		range = (0, 10),
		step = 0.1,
		decimal = 1,
		default = 0.5,
		display = "Vertical line width:"
	),
	RNASuiteParameter(
		key = 'show_hline',
		type = 'bool',
		default = True,
		display = "Show horizontal line:"
	),
	RNASuiteParameter(
		key = 'hline_type',
		type = 'list',
		options = ['dashed', 'dotted', 'dotdash', 'longdash', 'twodash', 'solid'],
		default = 'dashed',
		display = "Horizontal line type:"
	),
	RNASuiteParameter(
		key = 'hline_color',
		type = 'color',
		default = '#000000',
		display = "Horizontal line color:"
	),
	RNASuiteParameter(
		key = 'hline_width',
		type = 'float',
		range = (0, 10),
		step = 0.1,
		decimal = 1,
		default = 0.5,
		display = "Horizontal line width:"
	),
	RNASuiteParameter(
		key = 'theme_name',
		type = 'list',
		options = ['bw', 'classic', 'linedraw',
			'minimal', 'void', 'light', 'grey',
			'gray', 'dark'
		],
		display = 'Plot theme:'
	),
	RNASuiteParameter(
		key = 'base_size',
		type = 'int',
		range = (1, 100),
		step = 1,
		default = 11,
		display = 'Base size:'
	),
	RNASuiteParameter(
		key = 'legend_position',
		type = 'list',
		options = ['top', 'bottom', 'right'],
		default = 'top',
		display = 'Legend position:'
	),
	RNASuiteParameter(
		key = 'x_limit',
		type = 'limit',
		range = (-1000, 1000),
		step = 1,
		decimal = 2,
		default = (0, 0),
		display = 'X limit:'
	)
)



