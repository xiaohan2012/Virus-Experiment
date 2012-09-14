from urllib2 import urlopen
from codecs import open
from collections import defaultdict

from threading import Thread
from Queue import Queue

    
class DownloadThread(Thread):
    def __init__(self,queue):
        Thread.__init__(self)
        self.q = queue

    def _download(self,url,save_path,data={},headers={}):
        f = open(save_path,"w","utf8").write(urlopen(url).read())
        f.close()
        return True

    def run(self):
        while True:
            url,save_path= self.q.get()
            self._download(url,save_path)
            self.q.task_done()

class BatchDownloader():
    def __init__(self,url_path_list = [],thread_count = 10):
        self.thread_count = thread_count
        self.q = Queue()

        for url,path in url_path_list:
            self.q.put((url,path))

        self.is_started = False

    def addTasks(self,url_path_list):
        """task can be added halfway"""
        for url,path in url_path_list:
            self.q.put((url,path))
                    
    def start(self):
        """start downloading"""
        self.is_started = True
        self.thread_list = [DownloadThread(self.q) for i in xrange(self.thread_count)]
        for thread in self.thread_list:
            thread.setDaemon(True)
            thread.start()
        self.q.join()

    def isStarted(self):
        return self.is_started
    

                
class PageInformationFetcher(Thread):
    """fetch the download link contained in the webpage
    """
    def __init__(self,page_class , task_q , basket = [] , interested_info = {}):
        Thread.__init__(self)
        self.page_class = page_class
        self.basket = basket
        self.task_q = task_q
        self.intersted_info = interested_info

    def run(self):
        while True:
            page_url = self.task_q.get()

            page = self.page_class(page_url)
            print page_url
            page.prepare(page_url)

            for info_name,func_name in interested_info.items():
                info = getattr(page,func_name)()
                self.basket[info_name].append(info)
            #url  =page.getDownloadUrl()#fetch the url

            #self.url_basket.append(url)#put it in the basket

class BatchPageInformationFetcher(Thread):
    def __init__(self,page_class , page_urls=[],fetcher_count = 5,interested_info = {}):
        self.task_q = Queue()
        self.fetcher_count = fetcher_count
        for url in page_urls:
            self.task_q.put(url)

        self.page_class = page_class#the class that represents the page where interested information hides in
        
        self.basket = defaultdict(list)#container to hold the information

        self.is_start = False#indicate whether fetching has begun
        self.interested_info = interested_info #a dict hold, the interested info names and corresponding fetching method

    def isStarted(self):
        return self.is_start

    def addPageUrls(self,page_urls):
        for url in page_urls:
            self.task_q.put(url)
        
    def fetch(self):
        self.fetchers = [PageInformationFetcher(self.page_class , self.task_q , self.basket , self.interested_info)\
                             for i in xrange(self.fetcher_count)]
        
        for f in self.fetchers:
            f.setDaemon(True)
            f.start()
    def getAllInformation(self):
        self.task_q.join()
        return self.url_basket

if __name__ == "__main__":
    #prepare the information fetcher
    from pdbbank_crawler import SearchResultPage,ItemPage,load_page_content

    interested_info = {
        "download_url":"getDownloadUrl",
        "pdb_id":"getPDBId",
        "description":"getDescription",
        }
    info_fetcher = BatchPageInformationFetcher(ItemPage,interested_info = interested_info)

    url_template = "http://www.rcsb.org/pdb/results/results.do?gotopage=%d&qrid=4D35F462&tabtoshow=Current"

    srp = SearchResultPage(load_page_content(url_template %1))
    srp.prepare(url_template)

    page_urls = srp.getItemUrls()
    info_fetcher.addPageUrls( page_urls )
    info_fetcher.fetch()

    while srp.gotoNextPage():
        page_urls = srp.getItemUrls()[:2]
        info_fetcher.addPageUrls( page_urls )
        break
    info = info_fetcher.getAllInformation() 
    print info
