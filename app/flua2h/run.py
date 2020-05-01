from .tasks import flua2h
import json
if __name__ == '__main__':
	result = flua2h.delay('h5n1','HA','Chicken')
	data=result.get(timeout=5000)
	print(data)
	#print('Task result:'+str(data))
