import re
from pyquery import PyQuery as pq

from urllib2 import urlopen , Request
from urlparse import urlsplit, urljoin
from collections import defaultdict

from error import NoMoreNextError

def load_page_content(url):
    headers = {
            "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Charset":"ISO-8859-1,utf-8;q=0.7,*;q=0.3",
            "Accept-Language":"en-US,en;q=0.8",
            "Cache-Control":"max-age=0",
            "Connection":"keep-alive",
            "Cookie":"JSESSIONID=DB6E894842DAFECC4FF1899C8392D822; __utma=230211871.948318794.1347599384.1347599384.1347599384.1; __utmb=230211871.2.10.1347599384; __utmc=230211871; __utmz=230211871.1347599384.1.1.utmcsr=(direct)|utmccn=(direct)|utmcmd=(none)",
            "Host":"www.rcsb.org",
            "User-Agent":"Mozilla/5.0 (X11; Linux i686) AppleWebKit/537.1 (KHTML, like Gecko) Chrome/21.0.1180.79 Safari/537.1",
        }
    req = Request(url , headers = headers)     
    res = urlopen(req)
    return res.read()
    
class ItemPage(pq):
    def __init__(self,*args, **kwargs):
        pq.__init__(self , *args , **kwargs)#load the page
    
    def prepare(self,url):
        res = urlsplit(url)
        self.url = url
        self.baseurl = "%s://%s" %(res.scheme , res.netloc)

    def getPDBId(self):
        return re.findall(r'structureId=(\w+)',self.url)[0]

    def getDescription(self):
        key_col = self("#widget_se_unitMDAU .contentNB").find("td.mdauCategory")
        value_col = self("#widget_se_unitMDAU .contentNB").find("td.mdauData")
        desc = defaultdict(dict)
        for key , value in zip(key_col , value_col):
            key = pq(key).text().lower().strip(":")
            value = pq(value).text().strip()
            if key == "molecule":
                cur_molecule = value
            else:                
                desc[cur_molecule][key]=value
        return desc            
    
    def getDownloadUrl(self):
        try:
            return self._prepare_url(pq(self("#se_downloadFiles").find("a:contains('PDB File (Text)')")).attr("href"))
        except:
            return self._prepare_url(pq(self("#se_downloadFiles").find("a")[2]).attr("href"))

    def _prepare_url(self , url):
        """prepend the url with the baseurl"""
        return urljoin(self.baseurl , url)
class SearchResultPage(pq):
    def __init__(self,*args , **kwargs):
        pq.__init__(self,*args , **kwargs)
        self.page_index = 1

    def prepare(self,url_template):
        """
        set the url property using a separate function
        pquery is somewhat odd, args and kwargs should be used in the __init__ function
        """
        res = urlsplit(url_template)
        self.baseurl = "%s://%s" %(res.scheme , res.netloc)

        self.url_template = url_template 

        self._determine_last_page_index()

    def _determine_last_page_index(self):
        """anaylze the page and find/set the last page index"""
        tds = pq(self(".lb_head3")[-1]).find("td")
        last_url = pq(pq(tds[-1]).find("a")[-1]).attr("href")
        self.last_index = int(re.findall(r"gotopage=(\d+)" , last_url)[0])


    def _prepare_url(self , url):
        """prepend the url with the baseurl"""
        return urljoin(self.baseurl , url)

    def _get_next_url(self):
        if self.page_index + 1 <= self.last_index:
            self.page_index += 1
            return self.url_template %self.page_index
        else:
            return ''            

    def getItemUrls(self):
        return [self._prepare_url( pq(a).attr("href"))  for a in self.find("a.qrb_structid")]


    def gotoNextPage(self):
        url = self._get_next_url()
        if url:
            print url
            pq.__init__(self,load_page_content(url))
            return True
        else:
            return False            

                    
        
if __name__ == "__main__":
    from pprint import pprint
    url_template = "http://www.rcsb.org/pdb/results/results.do?gotopage=%d&qrid=4D35F462&tabtoshow=Current"

    srp = SearchResultPage(load_page_content(url_template %1))
    srp.prepare(url_template)

    urls = srp.getItemUrls()
    while srp.gotoNextPage():
        srp.getItemUrls()
    
    #pdb_url = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1BJ1"
    #item = ItemPage(pdb_url)
    #item.prepare(pdb_url)
    #print item.getPDBId()
    #print item.getDescription()
    #print item.getDownloadUrl()
