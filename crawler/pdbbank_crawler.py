from error import NoMoreNextError
from pyquery import PyQuery as pq

from urlparse import urlsplit, urljoin

class BankItemPage(pq):
    def __init__(self,*args , **kwargs):
        pq.__init__(self,*args , **kwargs)
    
    def getPDBId(self):
        return
    
    def getChainName(self):
        return
    
    def getDescription(self):
        return

class SearchResultPage(pq):
    def __init__(self,*args , **kwargs):
        url = args[0]
        res = urlsplit(url)
        self.baseurl = "%s://%s" %(res.scheme , res.netloc)
        print self.baseurl
        pq.__init__(self,*args , **kwargs)
    
    def getItemUrls(self):
        for a in self.find("a.qrb_structid"):
            print a
#print pq(a).text(), urljoin(self.baseurl , pq(a).attr("href"))
#yield  pq(a).text(), urljoin(self.baseurl , pq(a).attr("href"))

    def _get_next_url(self):
        return

    def _get_prev_url(self):
        return

    def gotoNextPage(self):
        url = self._get_next_url()
        if url:
            pq.__init__(self,url)
        else:
            raise NoMoreNextError
    
    def gotoPrevPage(self):        
        url = self._get_prev_url()
        if url:
            pq.__init__(self,url)
        else:
            raise NoMorePrevError
                
if __name__ == "__main__":
    search_result_url = "http://www.rcsb.org/pdb/results/results.do?qrid=10F3C09A&tabtoshow=Current"
    srp = SearchResultPage(search_result_url)
    from pprint import pprint
    pprint( srp.getItemUrls() )
            
