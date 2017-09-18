# ******************************************************************************
# pysimm.utils module
# ******************************************************************************
#
# ******************************************************************************
# License
# ******************************************************************************
# The MIT License (MIT)
#
# Copyright (c) 2016 Michael E. Fortunato, Coray M. Colina
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import print_function
from collections import Sequence
import itertools
from pprint import pprint

from pysimm import debug_print


class PysimmError(Exception):
    pass


class Container(object):
    """pysimm.utils.Container

    Abitrary container object that returns None if trying to access an attribute that does not exist
    """
    def __getattr__(self, name):
        return None


class ItemContainer(Sequence):
    """pysimm.utils.ItemContainer

    Container object intended to organize Item objects. Arbitrary attributes can be set using keyword arguments. Underlying data structure is a dictionary where the key is referred to as a tag, and the value should be an Item object. Item.tag should equal the key for the object in the dictionary.
    """
    def __init__(self, _dict=None, **kwargs):
        self._dict = _dict or {}
        self.count = len(_dict) if _dict else 0
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        for v in self._dict.values():
            yield v

    def __getitem__(self, slice_):
        if not isinstance(slice_, slice) and isinstance(slice_, int):
            if slice_ >= 0:
                return self._dict.get(slice_)
            else:
                return self._dict.get(self.count + slice_ + 1)
        elif isinstance(slice_, slice):
            data = []
            start, stop, step = slice_.indices(len(self._dict))
            for i in xrange(start+1, stop+1, step):
                item = self._dict.get(i)
                if item:
                    data.append(item)
            return data
        else:
            return None

    def add(self, _item):
        if _item.tag is not None and self._dict.get(_item.tag) is None:
            self._dict[_item.tag] = _item
            self.count += 1
        elif _item.tag is None and not self._dict.get(self.count+1):
            self.count += 1
            self._dict[self.count] = _item
            self._dict[self.count].tag = self.count
        else:
            debug_print('cannot add at index %s' % (self.count+1))
            return None
        return _item

    def get(self, *args, **kwargs):
        name = None
        tags = kwargs.get('tags')
        startswith = kwargs.get('startswith')
        first = kwargs.get('first')
        order = kwargs.get('order')
        item_wildcard = kwargs.get('item_wildcard') if kwargs.has_key('item_wildcard') else 'X'
        query_wildcard = kwargs.get('query_wildcard')
        improper_type = kwargs.get('improper_type')
        if len(args) == 1:
            name = args[0]
            if name == 'all':
                return self._dict.values()
        if (len(args) != 1 and kwargs.get('startswith') is None or
                kwargs.get('startswith') and kwargs.get('tags')):
            print('improper usage of system.ItemContainer.get')
            return []
        found = []
        for item in self:
            if name is not None and (item.name == name or
                                     (not order and item.rname == name)):
                found.append(item)
                if first:
                    return found

            elif name:
                if compare(name, item.name, query_wildcard=query_wildcard,
                           item_wildcard=item_wildcard,
                           improper_type=improper_type, order=order):
                    found.append(item)
                    if first:
                        return found
            elif startswith:
                if item.name and item.name.startswith(startswith):
                    found.append(item)
                elif (not item.name and item.type and item.type.name and
                        item.type.name.startswith(startswith)):
                    found.append(item)
                if first:
                    return found
            elif tags and isinstance(tags, list):
                if item.tag in tags:
                    found.append(item)
                    if first:
                        return found
        return sorted(found, key=lambda x: x.name.count('X'))

    def remove(self, index, update=True):
        if index == 'all':
            self._dict = {}
            self.count = 0
        if self._dict.get(index):
            if not update or index == self.count:
                del self._dict[index]
            else:
                for i in xrange(index, self.count):
                    assert self._dict[i+1] is not None
                    self._dict[i] = self._dict[i+1]
                    self._dict[i].tag = i
                del self._dict[self.count]
            self.count -= 1


class Item(object):
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if v == 'TRUE':
                v = True
            elif v == 'FALSE':
                v = False
            setattr(self, k, v)

    def __getattr__(self, name):
        return None

    def copy(self):
        i = type(self)()
        for k, v in vars(self).items():
            setattr(i, k, v)
        if i.tag:
            del i.tag
        return i

    def set(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


def compare(query, item, query_wildcard=None, item_wildcard='X', order=False,
            improper_type=False):
    query = query.split(',')
    item = item.split(',')

    if len(query) != len(item):
        print('cannot compare names with different lengths: %s %s'
              % (query, item))
        return False
    match = []
    if improper_type:
        if query[0] != item[0]:
            return False
        for q_perm in itertools.permutations(query[1:]):
            for i_perm in itertools.permutations(item[1:]):
                if compare(','.join(q_perm), ','.join(i_perm)):
                    return True
    else:
        for q, i in itertools.izip(query, item):
            if q == i:
                match.append(True)
            elif (item_wildcard and i == item_wildcard) or \
                    (query_wildcard and q == query_wildcard):
                match.append(True)
            else:
                match.append(False)
        if False not in match:
            return True

        if not order:
            match = []
            rquery = list(reversed(query))
            for q, i in itertools.izip(rquery, item):
                if q == i:
                    match.append(True)
                elif (item_wildcard and i == item_wildcard) or \
                        (query_wildcard and q == query_wildcard):
                    match.append(True)
                else:
                    match.append(False)
            if False not in match:
                return True

    return False
