from mpi4py import MPI
import numpy as np
import random
import math

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
stat = MPI.Status()

# set grid shape and number of subROWS, define origin, num generations
ROWS = 101
COLS = 101
if size > ROWS:
    print('Not enough ROWS')
    exit()
subROWS = ROWS // size + 2
mRows = size * (ROWS // size)
mCols = COLS
origin = (mRows//2, mCols//2)
generations = 50

# set model constants
k1, k2, k3, k4 = (0.74, 0.2, 0.4, 0.4)
phi = ROWS * COLS 
rho = 3.85 

# cell types: N, C, E, D (normal, cancerous, effector, dead)

# function definitions
def sumGrid(M, cType):
    # number of cells of certain type in entire grid
    return sum(sum(M==cType))

def distOrigin(p, orgn):
    # get euclidean distance from point (x,y) on grid to origin
    return math.sqrt((p[0] - orgn[0])**2 + (p[1] - orgn[1])**2)

def getQuadrant(coord, orgn):
    if coord[0] <= orgn[0] and coord[1] <= orgn[1]:
	return 'II'
    elif coord[0] <= orgn[0] and coord[1] > orgn[1]:
	return 'I'
    elif coord[0] > orgn[0] and coord[1] <= orgn[1]:
	return 'III'
    else:
	return 'IV'

def getCoords(r, c, sRows, mRows, rnk):
    # get cell position in mRows x mCols matrix 
    x = rnk * mRows + r - 1 
    return (x, c) 

def getTypeSumRnPrime(grid, subRows):
    # get num cancer cells, num abnormal cells, cancer cell total dist from origin
    c = sumGrid(grid[1:subRows-1,:],'C')
    e = sumGrid(grid[1:subRows-1,:],'E')
    d = sumGrid(grid[1:subRows-1,:],'D')
    nPrime = c + e + d
    R = 0
    for r in range(1,subRows-1):
	for c in range(0,COLS):
	    if grid[r,c] == 'C':
		coord = getCoords(r, c, subRows, mRows, rank)		
	 	R += distOrigin(coord, origin)
    return (c, nPrime, R)	

def densityDevelopment(nprime, r): 
    # calculate phi(t) based on num abnormal cells & total cancer cell dist from origin
    R = r / nprime
    return (nprime / (R**2))

def msgUp(subGrid):
    # send and receive rows with rank+1
    comm.send(subGrid[subROWS-2,:],dest=rank+1)
    subGrid[subROWS-1,:]=comm.recv(source=rank+1)
    return 0

def msgDn(subGrid):
    # send and receive rows with rank-1
    comm.send(subGrid[1,:],dest=rank-1)
    subGrid[0,:] = comm.recv(source=rank-1)
    return 0
 
def computeGrowth(subGrid, mProb, dense):
    # computes growth changes in each worker's subgrid for each generation
    newGrid = np.copy(subGrid)
    for r in range(0,subROWS):
	for c in range(1,COLS-1):    
	    if subGrid[r,c] == 'C':
		if random.random() < mProb:
		    newGrid = mitosis(r, c, newGrid, origin, dense) 
		elif random.random() < k2:
		    newGrid[r,c] = 'E'
	    elif subGrid[r,c] == 'E':
		if random.random() < k3:
		    newGrid[r,c] = 'D'
	    elif subGrid[r,c] == 'D':
		if random.random() < k4:
		    newGrid[r,c] = 'N'
    return newGrid

def mitosis(r, c, newGrid, origin, dense):
    # compute growth of cancerous cells through mitosis
    up = (r-1, c)
    rt = (r, c+1)
    dn = (r+1, c)
    lt = (r, c-1)
    abNormal = ('C','E','D')
    denseMap = {'I':[up,rt],'II':[lt,up],'III':[dn,lt],'IV':[rt,dn]}
    notDenseMap = {'I':[dn,lt],'II':[rt,dn],'III':[up,rt],'IV':[lt,up]}
    zeroDenseMap = {'I':[rt,rt],'II':[lt,lt],'III':[dn,lt],'IV':[rt,dn]}
    zeroNotDenseMap = {'I':[dn,lt],'II':[rt,dn],'III':[rt,rt],'IV':[lt,lt]}
    lastDenseMap = {'I':[up,rt],'II':[lt,up],'III':[lt,lt],'IV':[rt,rt]}
    lastNotDenseMap = {'I':[lt,lt],'II':[rt,rt],'III':[up,rt],'IV':[lt,up]} 
    choice = random.random()
    coord = getCoords(r, c, subROWS, mRows, rank) 
    quadrant = getQuadrant(coord, origin)
    if dense:
	if r == 0:
	    if choice < 0.5 and newGrid[zeroDenseMap[quadrant][0]] not in abNormal:
                newGrid[zeroDenseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[zeroDenseMap[quadrant][1]] not in abNormal:
                newGrid[zeroDenseMap[quadrant][1]] = 'C' 
	elif r == (subROWS-1):
	    if choice < 0.5 and newGrid[lastDenseMap[quadrant][0]] not in abNormal:
                newGrid[lastDenseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[lastDenseMap[quadrant][1]] not in abNormal:
                newGrid[lastDenseMap[quadrant][1]] = 'C'
	else:
            if choice < 0.5 and newGrid[denseMap[quadrant][0]] not in abNormal:
                newGrid[denseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[denseMap[quadrant][1]] not in abNormal:
                newGrid[denseMap[quadrant][1]] = 'C'
    else:
	if r == 0:
	    if choice < 0.5 and newGrid[zeroNotDenseMap[quadrant][0]] not in abNormal:
                newGrid[zeroNotDenseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[zeroNotDenseMap[quadrant][1]] not in abNormal:
                newGrid[zeroNotDenseMap[quadrant][1]] = 'C'	
	elif r == (subROWS-1):
	    if choice < 0.5 and newGrid[lastNotDenseMap[quadrant][0]] not in abNormal:
                newGrid[lastNotDenseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[lastNotDenseMap[quadrant][1]] not in abNormal:
                newGrid[lastNotDenseMap[quadrant][1]] = 'C'
	else:
            if choice < 0.5 and newGrid[notDenseMap[quadrant][0]] not in abNormal:
                newGrid[notDenseMap[quadrant][0]] = 'C'
            elif choice >= 0.5 and newGrid[notDenseMap[quadrant][1]] not in abNormal:
                newGrid[notDenseMap[quadrant][1]] = 'C'
    return newGrid

# initialize starting subgrids
subGrid = np.empty([subROWS, COLS], dtype=object)
subGrid.fill('N')

# initialize center of large grid with cancer cell(s)
if rank == (size//2):
    if size % 2 == 1:
	subGrid[subROWS//2, COLS//2-1] = 'C'
	subGrid[subROWS//2, COLS//2+1] = 'C'
    else:
	subGrid[1, COLS//2-1] = 'C'
	subGrid[1, COLS//2+1] = 'C'

# show the full initial grid
Grid = comm.gather(subGrid[1:subROWS-1,:],root=0)
if rank == 0:
    Grid = np.vstack(Grid)
    print('Initial Grid')
    print(Grid[:])
    print('')

# model tumor growth for g generations
for g in range(generations):
    diags = getTypeSumRnPrime(subGrid, subROWS) 
    dense = False
    numCancers = comm.gather(diags[0],root=0)
    nPrime = comm.gather(diags[1],root=0)
    R = comm.gather(diags[2],root=0)
    if rank == 0:
	numCancers = sum(numCancers)
	nPrime = sum(nPrime)
	R = sum(R)
	dense = densityDevelopment(nPrime, R) > rho
    cancers = comm.bcast(numCancers,root=0) 
    dense = comm.bcast(dense,root=0)
    mitosisProb = k1 * (1 - cancers / phi) 
    newGrid = computeGrowth(subGrid, mitosisProb, dense) 
    if rank == 0:
	msgUp(newGrid)
    elif rank == size-1:
	msgDn(newGrid)
    else:
	msgUp(newGrid)
	msgDn(newGrid)
    subGrid = np.copy(newGrid)

# show final grid and number of each cell type present
finalGrid = comm.gather(subGrid[1:subROWS-1,:],root=0)
if rank == 0:
    finalGrid = np.vstack(finalGrid)
    print('Final Grid')
    print(finalGrid[:])
    print('Normal: ',sumGrid(finalGrid,'N'),'Cancer: ',sumGrid(finalGrid,'C'),'Effector: ',
	  sumGrid(finalGrid,'E'),'Dead: ',sumGrid(finalGrid,'D')) 
