from random import sample

class Point(object):
    def __init__(self,code,dist_mat):
        self.dist_mat = dist_mat
        self.code = code

    def dist_to(self,p):
        if self == p:
            return 0

        try:
            return self.dist_mat[self.code][p.code]
        except KeyError
            return self.dist_mat[p.code][self.code]

    def __hash__(self):
        return self.code

    def __eq__(self , other):
        return self.code == other.code

    def __repr__(self):
        return self.code

class Cluster(set):
    def __init__(self,init_point,dist_mat):
        set.__init__(self,[init_point])
        self.center = init_point
        self.dist_mat = dist_mat

    def dist_to(self,point):
        try:
            return dist_mat[self.center][point]
        except KeyError:
            return dist_mat[point][self.center]

    def receive_point(self.point):
        self.add(point)

    def update_center(self):
        min_dist_avg = float('inf')
        temp_center = None
        for p1 in self:
            temp = 0
            for p2 in self:
                if p1 == p2:
                    continue
                else:
                    temp += p1.dist_to(p2)
            avg = temp / float(len(self))                    
            if  avg < min_dist_avg:
                min_dist_avg = avg
                temp_center = p1
        print "center updated to",p2
        self.center = min_dist_avg



def my_kmeans(sample_codes , dist_mat , k , iter_cnt =  20):
    init_centroids = sample(sample_codes , k)
    clusters = [Cluster(c,dist_mat) for c in init_centroids]
    points = [Point(code,dist_mat) for code in sample_codes]

    for i in xrange(iter_cnt):
        print "%d's interation,assigning points "
        for p in points:
            min_dist = float('inf')
            closest_cluster = None
            for c in clusters:
                dist = c.dist_to(p)
                if dist < min_dist:
                    min_dist = dist
                    closest_cluster = c
            closest_cluster.receive_point(p)
        print "updating centers"
        for c in clusters:
            c.update_center()

def init_dist():
    d_ = {}

