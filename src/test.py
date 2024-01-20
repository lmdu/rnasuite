import os
import time
import rchitect
from rchitect import console
from rchitect.callbacks import def_callback
from rchitect.interface import process_events, peek_event

os.environ["RCHITECT_RETICULATE_CONFIG"] = "0"
# enable signal handlers
os.environ["RCHITECT_REGISTER_SIGNAL_HANDLERS"] = "1"
args = ["radian", "--quiet", "--no-restore-history"]


rchitect.init(args=args, register_signal_handlers=True)
rchitect.rcall('library', 'reticulate')
rchitect.rcall('library', 'httpgd')
rchitect.rcall('hgd')
rchitect.rcall('hgd_browse')

rchitect.reval('x<-c(1,2,3,4,5,6); y<-c(4,5,6,7,8,9); plot(x,y)')
rchitect.rprint(rchitect.robject('abc'))

#time.sleep(2)

peek_event()
process_events()

while True:
	if peek_event():
		process_events()
		time.sleep(0.01)