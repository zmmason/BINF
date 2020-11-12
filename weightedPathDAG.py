# !/usr/bin/env python3
# DAGs - a generalization of the alignment problem
# Problem 18 - Find the Longest Path in a DAG
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE
import sys


class LongestPathDAG:
    """
    Find a longest path between two nodes in an edge-weighted DAG.
    Given:  An integer representing the source node of a graph, followed by an integer representing
            the sink node of the graph, followed by an edge-weighted graph
    Return: The length of a longest path in the graph, followed by a longest path.
    """

    def __init__(self, contents):
        """Constructor: holds data to be shared and manipulated by all functions"""
        self.contents = contents  # complete input graph
        self.adjDict = {}  # adjacency dict to hold all nodes and their possible sink nodes
        self.edgeWeight = {}  # dict to hold the weight of each source->sink pair
        self.startNode = int(self.contents[0])  # initial source node of the input graph
        self.endNode = int(self.contents[1])  # final sink node of the input graph

    def nodeList(self):
        """Initialize the list of nodes and their corresponding sinks/paths."""
        for edge in self.contents[2:]:  # disregarding the source and final sink node
            nodeData = list(edge.split('->'))
            if int(nodeData[0]) not in self.adjDict:  # adds node info to dict
                self.adjDict[int(nodeData[0])] = 0  # initialized the nodes for updating
            else:
                continue
        return

    def getEdges(self):
        """Add all possible sink nodes to each source node in the graph to update the Adjacency dictionary."""
        for i in range(len(self.adjDict)):
            endNodes = []  # list to hold sink nodes for each node
            for edge in self.contents[2:]:  # disregarding the source and final sink node
                nodeData = list(edge.split('->'))  # splits graph data for (source | sink:weight)
                edgeData = list(nodeData[1].split(':'))  # splits graph sink data for (sink | weight)
                if int(nodeData[0]) == i:
                    endNodes.append(int(edgeData[0]))
            self.adjDict[i] = endNodes  # ads the sink nodes to its corresponding source node
        self.adjDict[self.endNode] = []  # adding the final sink node as a node
        return

    def weightList(self):
        """Create a list of weights for each corresponding 'node->sink' pair."""
        for edge in self.contents[2:]:  # disregarding the source and final sink node
            split = edge.split(":")  # splits graph data for (source->sink | weight)
            self.edgeWeight[split[0]] = int(split[1])  # adds source->sink : weight to dict
        return

    def longestPath(self, v, visited=None, path=None):
        """Get all possible paths in the input graph starting at the source node."""
        if visited is None:  # initializes visited list for each new node with source node
            visited = [v]
        if path is None:  # initialized new path for each node starting at the source node
            path = [v]
        paths = []  # list to hold all paths
        for t in self.adjDict[v]:
            if t not in visited and t in self.adjDict:  # verification of valid sink nodes
                t_path = path + [t]  # growing the path
                paths.append(t_path)
            else:
                break
            paths.extend(self.longestPath(t, visited[:], t_path))  # loops the function until all nodes are visited
        return paths

    def getWeight(self, pathString):
        """Get the maximum weight and longest path from the list of relevant paths."""
        maxWeight = 0  # holds the max weighted path in the list of relevant paths
        bestPath = ''  # holds the longest path associated with the max weight
        for path in pathString:
            totalWeight = 0  # initializing the weight of each path weight
            for weight in self.edgeWeight:
                if weight in path:  # checks the created string path for source->sink pair in weight dictionary
                    totalWeight += self.edgeWeight[weight]  # counting the weights of each source->sink pair
            if totalWeight > maxWeight:  # checks to make sure the path weight to become the new max path weight
                maxWeight = totalWeight  # sets the new maxWeight
                bestPath = path  # sets the new best path
            else:
                continue
        return maxWeight, bestPath


def main():
    """
    Find a longest path between two nodes in an edge-weighted DAG.
    Given:  An integer representing the source node of a graph, followed by an integer representing
            the sink node of the graph, followed by an edge-weighted graph
    Return: The length of a longest path in the graph, followed by a longest path.
    """
    # contents = ['0', '4', '0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3']  # sample input
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.strip())
    myDAQ = LongestPathDAG(contents)
    myDAQ.nodeList()  # initialized the adjacency dict with each node
    myDAQ.getEdges()  # completes the adjacency dict
    myDAQ.weightList()  # creates the source->sink wight association dictionary
    paths = myDAQ.longestPath(int(contents[0]))  # gets all possible paths in the input graph starting at source node
    completePaths = []  # initializing list to hold paths that start at and stop at the specified nodes - source/sink
    for path in paths:
        if path[0] == int(contents[0]) and path[-1] == int(contents[1]):  # verifies source and sink node are present
            completePaths.append(path)
    stringCompletePaths = []  # initializing the list of relevant paths as string using '->' notation
    for path in completePaths:
        pathString = contents[0]  # adds the source node
        for node in path[1:]:
            pathString += '->' + str(node)  # adds each item in the path to the string with the proper format
        stringCompletePaths.append(pathString)
    maxWeight, bestPath = myDAQ.getWeight(stringCompletePaths)  # gets the max weight and the associated longest path
    print('{}\n{}'.format(maxWeight, bestPath))


if __name__ == "__main__":
    main()
