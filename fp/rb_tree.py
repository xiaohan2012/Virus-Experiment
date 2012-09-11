BLACK = 0
RED = 1

class tree_node(object):
    def __init__(self,  key = None,  data = None,  color = RED):
        self.key = key
        self.data = data
        self.parent = self.left_child = self.right_child = None
        self.color = color
        self.nonzero = 1
        self.count = 1

    def __nonzero__(self):
        return self.nonzero

    def is_red(self):
        if self.color == RED:
            return 1
        else:
            return 0

    def is_black(self):
        if self.color == BLACK:
            return 1
        else:
            return 0

    def append_data(self,  additional_data):
        for elem in self.data:
            if additional_data == elem:
                print "Item already in data list... Ignoring"
            else:
                data.append(additional_data)
                self.count += 1

    def delete_data(self,  del_data):
        if del_data in self.data:
            self.data.remove(del_data)
            self.count -= 1



class red_black_tree(object):
    """
    Class for building and modifying binary red black tree
    """
    def __init__(self):
        self.sentinel = tree_node()
        self.sentinel.left = self.sentinel.right = self.sentinel
        self.sentinel.color = BLACK
        self.sentinel.nonzero = 0
        self.root = self.sentinel
        self.elements = 0

    def compare_keys(self,  key_x,  key_y):
        return cmp(key_x,  key_y)

    def rotate_left(self, node_x):

        node_y = node_x.right

        node_x.right = node_y.left
        if node_y.left != self.sentinel:
            node_y.left.parent =node_x

        if node_y != self.sentinel:
            node_y.parent = node_x.parent
        if node_x.parent:
            if node_x == node_x.parent.left:
                node_x.parent.left = node_y
            else:
                node_x.parent.right = node_y
        else:
            self.root = node_y

        node_y.left = node_x
        if node_x != self.sentinel:
            node_x.parent = node_y


    def rotate_right(self, node_x):
        
        node_y = node_x.left

        node_x.left = node_y.right
        if node_y.right != self.sentinel:
            node_y.right.parent = node_x

        if node_y != self.sentinel:
            node_y.parent = node_x.parent
        if node_x.parent:
            if node_x == node_x.parent.right:
                node_x.parent.right = node_y
            else:
                node_x.parent.left =node_y
        else:
            self.root = node_y

        node_y.right = node_x
        if node_x != self.sentinel:
            node_x.parent = node_y


    def insert_rebalance(self, node_x):

        while node_x != self.root and node_x.parent.is_red():
            if node_x.parent == node_x.parent.parent.left:
                node_y = node_x.parent.parent.right

                if node_y.is_red():
                    node_x.parent.color = BLACK
                    node_y.color = BLACK
                    node_x.parent.parent.color = RED
                    node_x = node_x.parent.parent

                else:
                    if node_x == node_x.parent.right:
                        node_x =node_x.parent
                        self.rotate_left(node_x)

                    node_x.parent.color = BLACK
                    node_x.parent.parent.color = RED
                    self.rotate_right(node_x.parent.parent)
            else:
                node_y = node_x.parent.parent.left

                if node_y.is_red():
                    node_x.parent.color = BLACK
                    node_y.color = BLACK
                    node_x.parent.parent.color = RED
                    node_x = node_x.parent.parent

                else:
                    if node_x == node_x.parent.left:
                        node_x = node_x.parent
                        self.rotate_right(node_x)

                    node_x.parent.color = BLACK
                    node_x.parent.parent.color = RED
                    self.rotate_left(node_x.parent.parent)

        self.root.color = BLACK


    def insert_node(self, key, data):
        
        current = self.root
        parent = None
        
        while current != self.sentinel:
            rc = self.compare_keys(key,  current.key)
            if rc == 0:
                current.append_data(data)
                return current
            parent = current
            if rc < 0:
                current = current.left
            else:
                current = current.right

        new_node = tree_node(key, data)
        new_node.left = new_node.right = self.sentinel
        new_node.parent = parent

        self.elements = self.elements + 1

        if parent:
            if self.compare_keys(key, parent.key) < 0:
                parent.left = new_node
            else:
                parent.right = new_node
        else:
            self.root = new_node

        self.insert_rebalance(new_node)
        return new_node


    def delete_rebalance(self, node_x):
        
        while node_x != self.root and node_x.is_black():
            if node_x == node_x.parent.left:
                node_w = node_x.parent.right
                if node_w.is_red():
                    node_w.color = BLACK
                    node_x.parent.color = RED
                    self.rotate_left(node_x.parent)
                    node_w = node_x.parent.right

                if node_w.left.is_black() and node_w.right.is_black():
                    node_w.color = RED
                    node_x = node_x.parent
                else:
                    if node_w.right.is_black():
                        node_w.left.color = BLACK
                        node_w.color = RED
                        self.rotate_right(node_w)
                        node_w = node_x.parent.right

                    node_w.color = node_x.parent.color
                    node_x.parent.color = BLACK
                    node_w.right.color = BLACK
                    self.rotate_left(node_x.parent)
                    node_x = self.root

            else:
                node_w = node_x.parent.left
                if node_w.is_red():
                    node_w.color = BLACK
                    node_x.parent.color = RED
                    self.rotate_right(node_x.parent)
                    node_w = node_x.parent.left

                if node_w.right.is_black() and node_w.left.is_black():
                    node_w.color = RED
                    node_x = node_x.parent
                else:
                    if node_w.left.is_black():
                        node_w.right.color = BLACK
                        node_w.color = RED
                        self.rotate_left(node_w)
                        node_w = node_x.parent.left

                    node_w.color = node_x.parent.color
                    node_x.parent.color = BLACK
                    node_w.left.color = BLACK
                    self.rotate_right(node_x.parent)
                    node_x = self.root

        node_x.color = BLACK

    def delete_node(self,  convict,  data = None):
        
        if not convict or convict == self.sentinel:
            return
            
        if convict.count > 1 and not data: 
            convict.delete_data(data)
            return

        if convict.left == self.sentinel or convict.right == self.sentinel:
            y = convict
        else:
            y = convict.right
            while y.left != self.sentinel:
                y = y.left

        if y.left != self.sentinel:
            x = y.left
        else:
            x = y.right

        x.parent = y.parent
        if y.parent:
            if y == y.parent.left:
                y.parent.left = x
            else:
                y.parent.right = x
        else:
            self.root = x

        if y != convict:
            convict.key = y.key
            convict.data = y.data

        if y.color == BLACK:
            self.delete_rebalance(x)

        del y
        self.elements = self.elements - 1


    def delete_node_by_key(self,  key,  data = None):
        
        del_node = self.find_node(key)
        
        self.delete_node(del_node,  data)


    def find_node(self, key):

        current = self.root

        while current != self.sentinel:
            rc = self.compare_keys(key, current.key)
            if rc == 0:
                return current
            else:
                if rc < 0:
                    current = current.left
                else:
                    current = current.right

        return None


    def find_nodes_in_range(self,  down,  up):
        
        out = []
        current = self.root
        
        cmp_down = self.compare_keys(down, current.key)
        
        #Searching for lowermost item
        while current.left != self.sentinel and self.compare_keys(down, current.left.key) <= 0:
            if self.compare_keys(down, current.key) <= 0:
                current = current.left
            else:
                current = current.right
        #Trawersing right until reaching rightmost item within search criteria
        while current and current != self.sentinel and self.compare_keys(up, current.key) >= 0:
            out.append(current)
            current = self.next_node(current)
        return out



    def list_nodes(self):
        cur = self.first_node()
        result = []
        while cur:
            result.append(cur)
            cur = self.next_node(cur)
        return result


    def first_node(self):
        cur = self.root
        while cur.left:
            cur = cur.left
        return cur


    def last_node(self):
        cur = self.root
        while cur.right:
            cur = cur.right
        return cur


    def next_node(self, prev):
        cur = prev
        if cur.right:
            cur = prev.right
            while cur.left:
                cur = cur.left
            return cur
        while 1:
            cur = cur.parent
            if not cur:
                return None
            if self.compare_keys(cur.key, prev.key)>=0:
                return cur


    def prev_node(self, next):
        cur = next
        if cur.left:
            cur = next.left
            while cur.right:
                cur = cur.right
            return cur
        while 1:
            cur = cur.parent
            if cur is None:
                return None
            if self.__cmp(cur.key, next.key)<0:
                return cur
