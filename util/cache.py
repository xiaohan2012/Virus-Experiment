from collections import defaultdict

class BaseCache(dict):
    def set(self,key,value):
        raise NotImplementedError()

    def get(self,key):
        raise NotImplementedError()

        

