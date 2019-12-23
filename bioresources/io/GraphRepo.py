# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.utils.translation import gettext_lazy as __
from django.shortcuts import reverse

from neomodel import db
from bioresources.views import labelize
from collections import defaultdict


def title(obj):
    return obj["title"] if "title" in obj else obj["name"]


def _rid(obj, label):
    return label + "_" + str(obj["rid"] if "rid" in obj else obj.id)


def _label(node):
    return [x.strip() for x in node.labels if x != "Resource"][0]


class GraphRepo():

    @staticmethod
    def get_neighborhood_ids(rid_org, label_src, label_dst, level=2):
        if level == 2:
            query = (
                'MATCH (r:%s {rid: %s})-[e1]-(n1)-[e2]-(n2) WHERE n1:%s OR n2:%s  RETURN  n1,n2')
            query = query % (label_src, rid_org, label_dst, label_dst)
        elif level == 1:
            query = (
                'MATCH (r:%s {rid: %s})-[e1]-(n1) WHERE n1:%s  RETURN  n1')
            query = query % (label_src, rid_org, label_dst)
        else:
            raise Exception("Invalid level")
        results, meta = db.cypher_query(query, {})

        addedn = []

        related_resources = []

        for row in results:
            for node in row:
                rid = node["rid"]
                if (_label(node) == label_dst) and (rid_org != rid) and (rid not in addedn):
                    related_resources.append(rid)
                    addedn.append(rid)

        return related_resources

    @staticmethod
    def get_neighborhood(rid_org, label, level=2):
        if level == 2:
            # query = ('MATCH (r:%s {rid: {rid}})-[e1]-(n1)-[e2]-(n2) WHERE NOT n1:Country RETURN  r,e1,e2,n1,n2' % label)
            query = ('MATCH (r:%s {rid: {rid}})-[e1]-(n1)-[e2]-(n2:Resource)  RETURN  r,e1,e2,n1,n2' % label)
        elif level == 1:
            query = ('MATCH (r:%s {rid: {rid}})-[e1]-(n1) RETURN  r,e1,n1' % label)
        else:
            raise Exception("Level not allowed!")
        results, meta = db.cypher_query(query, {"rid": rid_org})
        nodes = []
        edges = []
        added = []
        addedn = []

        related_resources = defaultdict(lambda: [])
        rid_org = label + "_" + str(rid_org)
        for row in results:
            for node_or_edge in row:
                if "Relationship" in str(node_or_edge):
                    edge = node_or_edge
                    s, e = edge.start_node, edge.end_node
                    _from, _to = _rid(s, _label(s)), _rid(e, _label(e))
                    if set([_from, _to]) not in added:
                        edges.append({"from": _from, "to": _to})
                        added.append(set([_from, _to]))
                else:
                    node = node_or_edge
                    label = _label(node)
                    rid = _rid(node, label)
                    if rid not in addedn:
                        nodes.append({"id": rid, "label": labelize(title(node)), "ntype": label})
                        addedn.append(rid)
                        if (rid_org != rid)  and ("Resource" in node.labels):
                            related_resources[label].append(
                                {"url": reverse('bioresources:%s_view' % label.lower(), args=[node["rid"]]),
                                 "title": title(node)})
        related_resources = dict(related_resources)

        return {"nodes": nodes, "edges": edges}, related_resources
