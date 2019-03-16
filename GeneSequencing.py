#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	
# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []

		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):

				if(j < i):
					s = {}
				else:
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
					if not banded:
						if len(sequences[i]) > align_length:
							seq1 = sequences[i][0:align_length]
						else:
							seq1 = sequences[i]
						if len(sequences[j]) > align_length:
							seq2 = sequences[j][0:align_length]
						else:
							seq2 = sequences[j]
						if seq1 == seq2:
							score = len(seq2) * MATCH
							alignment1 = alignment2 = seq2
						else:
							score, alignment1, alignment2 = self.edit_dist(seq1, seq2)
					else: #banded
						if len(sequences[i]) > align_length:
							seq1 = sequences[i][0:align_length]
						else:
							seq1 = sequences[i]
						if len(sequences[j]) > align_length:
							seq2 = sequences[j][0:align_length]
						else:
							seq2 = sequences[j]
						if seq1 == seq2:
							score = len(seq2) * MATCH
							alignment1 = alignment2 = seq2
						else:
							score, alignment1, alignment2 = self.banded_edit_dist(seq1, seq2)
###################################################################################################					
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()	
				jresults.append(s)
			results.append(jresults)
		return results

	def edit_dist(self, seq1, seq2):
		matrix = []
		for i in range(len(seq1) + 1):
			matrix.append([])
			for j in range(len(seq2) + 1):
				if i == 0 and j == 0:
					node = self.Node()
					node.backtrace = float('inf')
					node.cost = 0
					matrix[i].append(node)
				elif i == 0:
					node = self.Node()
					node.backtrace = (i, (j - 1))
					node.cost = matrix[i][j - 1].cost + INDEL
					matrix[i].append(node)
				elif j == 0:
					node = self.Node()
					node.backtrace = ((i - 1), j)
					node.cost = matrix[i - 1][j].cost + INDEL
					matrix[i].append(node)
				else:
					node = self.Node()
					node.set_node(seq1[i - 1], seq2[j - 1], i, j, matrix)
					matrix[i].append(node)

		score = matrix[len(matrix) - 1][len(matrix[1]) - 1].cost
		node = matrix[len(matrix) - 1][len(matrix[1]) - 1]
		while node.backtrace != float('inf'):
			if seq1[node.backtrace[0]] == seq2[node.backtrace[1]] or (node.backtrace[0] == (node.curr_node[0] - 1) and node.backtrace[1] == (node.curr_trace[1] - 1)): #match or substitution
				alignment1 = seq1[node.curr_trace[0]] + alignment1
				alignment2 = seq2[node.curr_trace[1]] + alignment2
			elif node.backtrace[0] == (node.curr_trace[0] - 1):
				alignment2 = "-" + alignment2
				alignment2 = seq2[node.curr_trace[1]] + alignment2
			elif node.backtrace[1] == (node.curr_trace[1] - 1):
				alignment1 = seq1[node.curr_trace[0]] + alignment1
				alignment1 = "-" + alignment1
		return score, alignment1, alignment2

	def banded_edit_dist(self, seq1, seq2):
		matrix = []
		for i in range(len(seq1) + 1):
			matrix.append([])
			for j in range(len(seq2) + 1):
				if i == 0 and j == 0:
					node = self.Node()
					node.backtrace = float('inf')
					node.cost = 0
					matrix[i].append(node)
				elif i == 0:
					node = self.Node()
					node.backtrace = (i, (j - 1))
					node.cost = matrix[i][j - 1].cost + INDEL
					matrix[i].append(node)
				elif j == 0:
					node = self.Node()
					node.backtrace = ((i - 1), j)
					node.cost = matrix[i - 1][j].cost + INDEL
					matrix[i].append(node)
				else:
					node = self.Node()
					node.set_node(seq1[i - 1], seq2[j - 1], i, j, matrix)
					matrix[i].append(node)

		score = matrix[len(matrix) - 1][len(matrix[1]) - 1].cost
		alignment1 = ""
		alignment2 = ""
		return score, alignment1, alignment2

	class Node:
		backtrace = float('inf')
		curr_trace = (0,0)

		def __init__(self):
			pass

		def set_node(self, char1, char2, x, y, matrix):
			self.curr_trace = (x, y)
			if char1 == char2:
				self.backtrace = ((x - 1), (y - 1))
				self.cost = matrix[self.backtrace[0]][self.backtrace[1]].cost + MATCH
			else:
				left = matrix[(x - 1)][ y].cost + INDEL
				up = matrix[x][(y - 1)].cost + INDEL
				diag = matrix[(x - 1)][(y - 1)].cost + SUB

				if left <= up and left <= diag:
					self.backtrace = ((x - 1), y)
					self.cost = left
				if up <= left and up <= diag:
					self.backtrace = (x, (y - 1))
					self.cost = up
				if  diag <= left and diag <= up:
					self.backtrace = ((x - 1), (y - 1))
					self.cost = diag

		def banded_set_node(self, char1, char2, x, y, matrix):
			if char1 == char2:
				self.backtrace = ((x - 1), (y - 1))
				self.cost = matrix[self.backtrace[0]][self.backtrace[1]].cost + MATCH
			else:
				left = matrix[(x - 1)][ y].cost + INDEL
				up = matrix[x][(y - 1)].cost + INDEL
				diag = matrix[(x - 1)][(y - 1)].cost + SUB

				if left <= up and left <= diag:
					self.backtrace = ((x - 1), y)
					self.cost = left
				if up <= left and up <= diag:
					self.backtrace = (x, (y - 1))
					self.cost = up
				if  diag <= left and diag <= up:
					self.backtrace = ((x - 1), (y - 1))
					self.cost = diag


