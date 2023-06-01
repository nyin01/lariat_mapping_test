from logging import Logger

def setup_logger(filename: str = None, log_dir=None) -> Logger:
	'''
	Creates and returns a Logger obj for logging in scripts with a specific setup
	If a filename is passed in, the logs will be written to that file. 
	If a filename isn't passed in, the logs will be directed to stdout
	'''
	from logging import basicConfig, getLogRecordFactory, getLogger, setLogRecordFactory
	from time import asctime, localtime, strftime
	from os.path import join

	def formatted_record(*args, **kwargs):
		
		record = base_factory(*args, **kwargs)
		record.msg = str(record.msg)

		current_time = strftime('%I:%M:%S %m/%d')

		mili_elapsed = record.relativeCreated	#The amount of time since logging was imported that dummy was created in miliseconds
		secs_elapsed = mili_elapsed/1000
		time_elapsed = f'{int(secs_elapsed // 3600)}h:{int(secs_elapsed // 60 % 60)}m:{int(secs_elapsed % 60)}s'
		
		timestamp = current_time + '  -  ' + time_elapsed
		timestamp += ' '*(31 - len(timestamp)) + '|  '
		record.msg = timestamp + record.msg
		return record

	if log_dir is not None:
		filename = join(log_dir, filename)
	
	base_factory = getLogRecordFactory()
	current_time = strftime('%m/%d/%y %I:%M:%S')

	setLogRecordFactory(formatted_record)
	# If filename was specified, write logs to the file
	if filename is not None:
		basicConfig(filename=filename,
					filemode='w',
					level=20,				#level "info" has value 20
					format='%(message)s')
	# If it wasn't, print logs in stdout
	else:
		basicConfig(level=20,				#level "info" has value 20
					format='%(message)s')

		
	log = getLogger()
	log.info(f'START TIME: {current_time}')
	return log