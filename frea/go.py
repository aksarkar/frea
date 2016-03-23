"""Recover most specific Gene Ontology terms from pathway enrichments.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import gzip
import sys

def parse_obo(file):
    """Parse Gene Ontology in Open Biomedical Ontology format.

    Based on
    https://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/

    """
    singletons = set(['id', 'name', 'namespace', 'def'])

    def cleanup_term(term):
        return {k: term[k][0] if k in singletons else term[k] for k in term}

    term = None
    for line in file:
        line = line.strip()
        if not line:
            continue
        elif line == '[Typedef]':
            term = None
        elif line == '[Term]':
            if term:
                yield cleanup_term(term)
            term = collections.defaultdict(list)
        elif term is not None:
            k, v = line.split(': ', maxsplit=1)
            real_v, *comment = v.split(' ! ', maxsplit=1)
            term[k].append(real_v)
    if term:
        yield cleanup_term(term)

def parse_enrichments(file):
    """Parse pathway enrichments"""
    for line in file:
        yield line.strip().split('\t')

def build_graph(file):
    """Return a dict of nodes and a dict of edges for Gene Ontology terms

    Nodes map term names to IDs. Edges follow is_a and other relationships.

    """
    nodes = {}
    edges = collections.defaultdict(list)
    for term in parse_obo(file):
        nodes[term['name']] = term['id']
        for parent in term.get('is_a', []):
            edges[term['id']].append(parent)
        for other in term.get('relationship', []):
            rel, parent = other.split()
            edges[term['id']].append(parent)
    return nodes, edges

def traverse(pathways, nodes, edges):
    """DFS coroutine

    Yields nodes when visited and parents when queued.

    """
    visited = set()
    for p in pathways:
        start = nodes[p]
        if start in visited:
            continue
        yield 'visit', start
        queue = [start]
        while queue:
            p = queue.pop()
            visited.add(p)
            for parent in edges[p]:
                if parent not in visited:
                    queue.append(parent)
                yield 'parent', parent

def resolve_enrichments(pathways, nodes, edges):
    """Find all enriched pathways for which no child is enriched"""
    fringe = set()
    subgraph = {k: pathways[k] for k in pathways if k in nodes}
    search = traverse(subgraph, nodes, edges)
    while True:
        try:
            op, node = search.send(None)
            if op == 'visit':
                fringe.add(node)
            elif op == 'parent':
                if node in fringe:
                    fringe.remove(node)
            else:
                raise ValueError
        except StopIteration:
            break
    lookup = {nodes[k]: k for k in nodes} 
    for p in fringe:
        yield lookup[p], pathways[lookup[p]]

if __name__ == '__main__':
    with gzip.open(sys.argv[1], 'rt') as f:
        nodes, edges = build_graph(f)

    with open(sys.argv[2]) as f:
        pathways = {k: v for k, v in parse_enrichments(f)}

    for pathway, genes in resolve_enrichments(pathways, nodes, edges):
        print(pathway, genes, sep='\t')
