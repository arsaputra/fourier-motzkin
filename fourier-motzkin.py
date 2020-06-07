from copy import deepcopy
from math import ceil
from math import floor
import sys

#Help routines as follows

def polyhedron_constructor(input_file):
    '''
    Reads the given file and constructs a Polyhedron Ax>=b

    @param  str  input_file        The name of the input file, along with filetype
    @return list [matrix,vector]   Matrix A and Vector b put together in a List
    '''
    file=open(input_file,'r')
    lines=file.readlines()
    file.close()
    content=[]
    for i in lines:
        j=i.replace('\n','')
        content.append(j)

    #Construct Matrix
        
    matrix=[]
    for i in content:
        if i[0]=='#' or i[0]=='A': 
            pass
        elif i[0]=='b':
            k=content.index(i[0])
            break
        else:
            z=[]
            split=i.split()
            for x in split:
                z.append(float(x))
            matrix.append(z)
    
    #Construct Vector
            
    v=content[k+1].split()
    vector=[[float(i)] for i in v]
    
    return [matrix,vector]


def matrix_constructor(input_file):
    '''
    Reads the given file and constructs a matrix

    @param  str           input_file  The name of the input file along with filetype
    @return list(list)    matrix      The matrix as list of list
    '''
    file=open(input_file,'r')
    lines=file.readlines()
    file.close()
    content=[]
    for i in lines:
        j=i.replace('\n','')
        content.append(j)
    matrix=[]
    for i in content:
        if i[0]=='#': 
            pass
        else:
            z=[]
            split=i.split()
            for x in split:
                z.append(float(x))
            matrix.append(z)
    return matrix



def projection(A,b):
    '''
    Gives the n-1 th Projection of P(A,b) with n=#Columns of A. We define this projection as P(A_new,b_new)

    @param   list(list)   A                               The Matrix A
    @param   list(list)   b                               The Vector b
    @return  list(list)   (A_new,b_new)  P(A_new,b_new) in a list, with duplicates eliminated
    '''
    #Rearranging Matrix
    #Combine Matrix and Vector so we have Extended Coeff. Matrix
    
    combined=[A[i]+b[i] for i in range(len(A))]
    
    #Sort ECM by last variable in Positive,Negative,Zero parts
    
    positive=[row for row in combined if row[-2]>0]
    negative=[row for row in combined if row[-2]<0]
    zero=[row for row in combined if row[-2]==0]
    combined=positive+negative+zero

    #Split back to original

    b=[[combined[i].pop()]for i in range(len(combined))]
    A=combined

    #We have now grouped A and b 
    
    n=len(A[0])

    #Positive,Negative,Zero
    
    if len(positive)!=0 and len(negative)!=0 and len(zero)!=0:
        m1=len(positive)-1           #PYTHON INDEX NUMBERS. m strich
        m2=m1+len(negative)          #m doppelstrich
        m3=len(A)-1                  #m

        A_new=[]
        b_new=[]
        
        for i in range(m1+1):
            for k in range(m1+1,m2+1):
                a2=[A[k][j]*A[i][-1]-A[i][j]*A[k][-1] for j in range(n-1)]
                b2=[b[k][0]*A[i][-1]-b[i][0]*A[k][-1]]
                A_new.append(a2)
                b_new.append(b2)

        for i in range(m2+1,m3+1):
            a2=[A[i][j] for j in range(n-1)]
            b2=b[i]
            A_new.append(a2)
            b_new.append(b2)

        return [A_new,b_new] 

    #Positive
    
    if len(positive)!=0 and len(negative)==0 and len(zero)==0:
        A_new=[[0 for i in range(len(A[0])-1)] for j in range(len(A))]
        b_new=[[0]for i in range(len(b))]
        return [A_new,b_new]

    #Negative
    
    if len(positive)==0 and len(negative)!=0 and len(zero)==0:
        A_new=[[0 for i in range(len(A[0]))] for j in range(len(A))]
        b_new=[[0]for i in range(len(b))]
        return [A_new,b_new]
    
    #Zero -> We're actually 'done'!
    
    if len(positive)==0 and len(negative)==0 and len(zero)!=0:
        for i in range(len(A)):
            A[i].pop()
        return [A,b]
                
    #Positive, Negative
    
    if len(positive)!=0 and len(negative)!=0 and len(zero)==0:
        m1=len(positive)-1           
        m2=m1+len(negative)                       

        A_new=[]
        b_new=[]
        
        for i in range(m1+1):
            for k in range(m1+1,m2+1):
                a2=[A[k][j]*A[i][-1]-A[i][j]*A[k][-1] for j in range(n-1)]
                b2=[b[k][0]*A[i][-1]-b[i][0]*A[k][-1]]
                A_new.append(a2)
                b_new.append(b2)

        return [A_new,b_new] 
        
    #Negative, Zero
    
    if len(positive)==0 and len(negative)!=0 and len(zero)!=0:
        m1=len(negative)
        m2=len(zero)+m1
        A_new=[]
        b_new=[]
        for i in range(m1,m2):
            a2=[A[i][j] for j in range(n-1)]
            b2=b[i]
            A_new.append(a2)
            b_new.append(b2)
        return [A_new,b_new]
    
    #Positive, Zero

    if len(positive)!=0 and len(negative)==0 and len(zero)!=0:
        m1=len(positive)
        m2=len(zero)+m1
        A_new=[]
        b_new=[]
        for i in range(m1,m2):
            a2=[A[i][j] for j in range(n-1)]
            b2=b[i]
            A_new.append(a2)
            b_new.append(b2)
        return [A_new,b_new]

def F_M(A,b,k):
    '''
    Does the Fourier-Motzkin Elimination k times

    @param   list(list)  A      The Matrix
    @param   list(list)  b      The Vector
    @param   int         k      The number of times we want to do the Fourier-Motzkin Elimination
    @return  list(list   [A,b]  The resulting Polyhedron in a list
    '''
    n=len(A[0])
    while n>k:
        n-=1
        A,b=projection(A,b)
    return [A,b]

def identity(n):
    '''
    Creates a nxn identity matrix
    
    @param   int         n  Dimension of the desired Matrix (nxn)
    @return  list(list)  I  A nxn Identity Matrix
    '''
    I=[[0 for j in range(n)]for i in range(n)]
    for i in range(n):
        I[i][i]=1
    return I

def minus(A):
    '''
    Outputs -A for a given matrix A

    @param   list(list)  A  The Matrix
    @return  list(list)  A  The Matrix -A
    '''
    for i in range(len(A)):
        for j in range(len(A[0])):
            A[i][j]=-A[i][j]
    return A

def transpose(A):
    '''
    Outputs the Transpose of a given Matrix A

    @param   list(list)  A   The Matrix
    @return  list(list)  AT  The transposed Matrix
    '''
    m,n=len(A),len(A[0])
    AT=[[0 for i in range(m)]for j in range(n)]
    for i in range(m):
        for j in range(n):
            AT[j][i]=A[i][j]
    return AT


def solve(LGS):
    '''
    Solving linear equation of the following form:
    a1*x=b1
    a2*x=b2
    .
    .
    an*x=bn

    @param   list(list)  LGS                      A Polyhedron depicting list containing a Matrix an a Vector
    @return  bool        False                    Output if given Inequality has no solution  
    @return  int         ceil(maxpos)             A feasible solution
    @return  int         floor(minneg)            A feasible solution
    @return  int         ((minneg+maxpos)/2)      A feasible solution
    '''
    A=LGS[0]
    b=LGS[1]
    if checkzero(A)==True and notpositive(b)==True:
        return 0
    else:
        combined=[A[i]+b[i] for i in range(len(A))]
        pos=[i for i in combined if i[-2]>0]
        neg=[i for i in combined if i[-2]<0]
        pos_answer=[i.pop() for i in pos]
        neg_answer=[i.pop() for i in neg]
        pos_inq=[pos_answer[i]/pos[i][0] for i in range(len(pos))]
        neg_inq=[neg_answer[i]/neg[i][0] for i in range(len(neg))]
        try:
            maxpos=max(pos_inq)
        except ValueError:
            maxpos=-float('inf')
        try:
            minneg=min(neg_inq)
        except ValueError:
            minneg=float('inf')
        if minneg<maxpos:
            return False
        else:
            if minneg==float('inf'):
                return ceil(maxpos)
            if maxpos==-float('inf'):
                return floor(minneg)
            if minneg!=float('inf') and maxpos!=-float('inf'):
                return (minneg+maxpos)/2


#CHECK ABOVE

            
def checkzero(A):
    for i in A:
        for j in i:
            if j!=0:
                return False
    return True

def notpositive(A):
    for i in A:
        for j in i:
            if j>0:
                return False
    return True

def checkboolvect(v):
    for i in v:
        if type(i)==bool:
            return True
    return False
    
def move(LGS):
    '''
    Suppose a linear inequality of the following form is given:
    a11+a12+...+a1n*x1=b1
    .
    .
    an1+an2+...+ann*xn=bn
    the routine isolates the variables x1,...,xn in the left hand side of the equation

    @param  list(list)  LGS                                                                               The linear inequality
    @return list(list)  [[coeff[i] for i in range(len(coeff))]]+[[vector[i] for i in range(len(coeff))]]  The linear inequality after moving the constants
    '''
    matrix=LGS[0]
    vector=LGS[1]
    coeff=[[i.pop()]for i in matrix]
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            vector[i][0]-=matrix[i][j]
    return [[coeff[i] for i in range(len(coeff))]]+[[vector[i] for i in range(len(coeff))]]


def compute_x(A,b):
    '''
    Outputs a vector x that satisfies Ax>=b
    We obtain this by taking the 1st projection of P(A,b) and finding the value of x1, substituting it to 2nd projection of P(A,b), and so on

    @param   list(list)  A  The Matrix
    @param   list(list)  b  The Vector
    @return  int         x  A feasible Solution to Ax>=b
    '''
    all_projections=[[A,b]]
    x=[[0]for i in range(len(A[0]))]
    for i in range(len(A[0])-1):
        A,b=projection(A,b)
        all_projections.append([A,b])
    all_projections.reverse()
    x[0][0]=solve(move(all_projections[0]))
    #print(x[0][0])
    for i in range(1,len(x)):
        for j in range(i):
            for k in range(len(all_projections[i][0])):
                all_projections[i][0][k][j]*=x[j][0]
        x[i][0]=solve(move(all_projections[i]))
        #print(x[i][0])
    return x

def create_unit_row(r,i):
    '''
     create unit row
     @param two numbers: r lenght of the row and i die position of 1
        @return list with a unit row
    '''
    unit_row=[]
    for t in range (r):
        unit_row.append(0)
    unit_row[i]=1
    return unit_row

def mul_number(b,a):
    '''
    multiplication of a vector with number
    @param a vector and a number
        @return a vector a*b
    '''
    b_mul=[]
    for c in b:
        b_mul.append(c*a)
    return b_mul

def mat_vect(A,b):
    C=[]
    for i in range(len(A)):
        C.append(0)
    for i in range(len(A)):
        for j in range(len(A[0])):
            C[i]+=A[i][j]*b[j]
    return C

def vector_minus(u,v):
    '''
    subtraction of two vectors
    @param two vectors u,v
        @return vector u-v
    '''
    minu=[]
    for i in range(len(u)):
        minu.append(u[i]-v[i])
    return minu

def build_M(A,b,k):

    '''
    @param A matrix
    @param b vector
    @param k which projection do we want to do
    @return [D,d,M] list of matrices and vectors
    finds the matrices M so that D=M*A and d=M*b (representing the projection)
    outputs a list containing the matrix D, vector d, and matrix M
    '''
    
    positive=[A.index(row) for row in A if row[k]>0]
    negative=[A.index(row) for row in A if row[k]<0]
    zero=[A.index(row) for row in A if row[k]==0]
    total=len(positive)*len(negative)+len(zero)
    M=null(total,len(A))

    #print(positive)
    #print('\n')
    #print(negative)
    #print('\n')
    #print(zero)
    #print('\n')
    
    R_1=[(s,t) for s in negative for t in positive]
    R_2=[(s,t) for s in negative for t in positive]
    R_1+=zero

    #print(R_1)
    #print('\n')
    #print(R_2)

    if len(R_1)!=0:
        for i in range(len(M)): 
                if R_1[i] in zero:
                    M[i][R_1[i]]=1
                if R_1[i] in R_2:
                    pos=R_1.index(R_1[i])
                    M[i][R_1[pos][0]]=A[R_1[pos][1]][k]
                    M[i][R_1[pos][1]]=-A[R_1[pos][0]][k]
        D,d=multiply(M,A),mat_vect(M,b)

    else:
        D=null(1,len(A[0]))
        d=[0]
        M=null(1,len(A))

    return [D,d,M]


def make_columns(B):
    '''
    rewrite a matrix B
    @param the matrix B (list of rows)
        @return the matrix B (list of columns)
    '''
    N=[]
    for i in range(len(B[0])):
        I=[]
        for m in range(len(B)):
            I.append(B[m][i])
        N.append(I)
    return N


def multiply(A,B):
        '''
        multiplication of two matrices
        @param matrices A,B
            @return product A*B
        '''
        B=make_columns(B)
        zm=(len(A))
        sm=(len(B))
        sma=(len(A[0]))
        S=[[] for i in range(len(A))]
        for o in range(zm):
            for g in range(sm):
                D=[]
                for y in range(sma):
                    D.append(A[o][y]*B[g][y])
                S[o].append(sum(D))
        return S

def mul_list(L):
    '''
    multiplication of all elemensts of L
    @param list of matrices(list of rows)
    @return produkt of all matrices in L
    '''
    J=[]
    if len(L)>1:
        J.append(multiply(L[0],L[1]))
        L.remove(L[1])
        L.remove(L[0])
        J=J+L
        return mul_list(J)
    else:
        return L


def einheitsvektor(n,i):
    '''
    Creates a unit vector, with 1 on the i-th Entry
    '''
    E=[[0] for i in range(n)]
    E[i][0]=1
    return E

def null(m,n):
    '''
    constructs a matrix of dimesion mxn. fills all entry initially with 0
    '''
    return [[0 for i in range(n)] for i in range(m)]
    
        
#Routines

def project(input_file,k,output_file):
    '''
    Outputs the Matrix and the vector written on the file of the desired name which depict the kth projection of P(A,b)
    (Details see attached file)

    @param   str  input_file   The name of the input file with extension
    @param   int  k            The number of times we want to apply FME
    @param   str  output_file  The name of the produced file with extension
    '''
    LGS=polyhedron_constructor(input_file) #Matrix,Vector
    A=LGS[0]
    b=LGS[1]
    A,b=F_M(A,b,k)
    f=open(output_file,'w+')
    f.write('# '+str(k)+'th Projection of Polyhedron P #')
    f.write('\n')
    f.write('A\n')
    for row in A:
        for i in row:
            f.write(str(i)+' ')
        f.write('\n')
    f.write('b\n')
    for row in b:
        f.write(str(row[0])+' ')
    f.close()

    
def image(input_polyhedron_file,input_matrix_file,output_file):
    '''
    Outputs the Polyhedron M*P
    (Details see attached file)

    @param  str  input_polyhedron_file  Name of the Polyhedron file with extension
    @param  str  input_matrix_file      Name of the Matrix file with extension
    @param  str  output_file            Name of the produced file with extension
    '''
    LGS=polyhedron_constructor(input_polyhedron_file)
    A=LGS[0] #Dimension m x n
    b=LGS[1]
    m=len(A)
    n1=len(A[0])

    M=matrix_constructor(input_matrix_file) #Dimension k x n

    k=len(M)
    n2=len(M[0])

    #check if the dimensions are compatible
    
    if n1!=n2:
        print('Error')

    #Construct giant matrix. Prepare components 
    
    else:
        I1=identity(k)
        I2=deepcopy(I1)
        I2=minus(I2)
        M2=deepcopy(M)
        M2=minus(M2)
        zero=[[0 for i in range(k)]for j in range(m)]

        BM1=[I2[i]+M[i]for i in range(k)]
        BM2=[I1[i]+M2[i] for i in range(k)]
        BM3=[zero[i]+A[i] for i in range(m)]
        GM=BM1+BM2+BM3 

        #Construct giant vector
        
        GV=[[0]for i in range(2*k)]+b
    
        N=k+n1
        #kth Projection of Dimension k+n
        
        while N>k:
            N-=1
            GM,GV=projection(GM,GV)
    
    f=open(output_file,'w+')
    f.write('# IMAGE OF POLYHEDRON P UNDER M #')
    f.write('\n')
    f.write('A\n')
    for row in GM:
        for i in row:
            f.write(str(i)+' ')
        f.write('\n')
    f.write('b\n')
    for row in GV:
        f.write(str(row[0])+' ')
    f.close()

def H_representation(input_file,output_file):
    '''
    Outputs a Polyhedron (represented by a Matrix and a Vector) that depicts conv{x1,...,xl}

    @param  str  input_file   The name of the input file with extension
    @param  str  output_file  The name of the output file with extension
    '''

    #Input is a file consisting of k lines representing xk, each having n components. We must estract the matrix and transpose it to have a nxk matrix
    
    M=matrix_constructor(input_file)
    M=transpose(M)
    n=len(M)
    k=len(M[0])
    
    #Dimension of M: n x k
    
    #Do analogue to image

    A=[]
    one=[1 for i in range(k)]
    minuss=[-1 for i in range(k)]
    ident=identity(k)
    A.append(one)
    A.append(minuss)
    A+=ident


    #Dimension of A: k+2 x k
    
    b=[[1],[-1]]
    zero=[[0] for i in range(k)]
    b+=zero


    #Dimension of b: k+2 x 1
    
    I1=identity(n)
    I2=deepcopy(I1)
    I2=minus(I2)
    M2=deepcopy(M)
    M2=minus(M2)
    zero=[[0 for i in range(n)]for j in range(k+2)]

    BM1=[I2[i]+M[i]for i in range(n)]
    BM2=[I1[i]+M2[i] for i in range(n)]
    BM3=[zero[i]+A[i] for i in range(k+2)]
    GM=BM1+BM2+BM3

    GV=[[0]for i in range(2*n)]+b
     
    N=n+k
        
    while N>n:
        N-=1
        GM,GV=projection(GM,GV)

    f=open(output_file,'w+')
    f.write('# CONVEX HULL OF x1,...,Xk #')
    f.write('\n')
    f.write('A\n')
    for row in GM:
        for i in row:
            f.write(str(i)+' ')
        f.write('\n')
    f.write('b\n')
    for row in GV:
        f.write(str(row[0])+' ')
    f.close()

def compute_x_or_y(input_file):
    '''
    Prints a vector x in P if P is nonempty and returns True, or prints a vector y with transpose(y)*A=0 and transpose(y)*b>0 if P is empty and returns False

    @param   str   input_file  The input file name with extension
    @return  bool  True        if P is nonempty
    @return  bool  False       if P is empty
    '''
    A,b=polyhedron_constructor(input_file)
    C,d=F_M(A,b,0)
    if notpositive(d)==True:
        
        #there exists an answer
        
        res=compute_x(A,b)
        new=[i[0] for i in res]
        print((True,new))
    else:
        #there is no answer. find infeasibility certificate
        
        spalte=len(A[0])
        A1,b1=deepcopy(A),deepcopy(b)
        b2=[i[0] for i in b1]
        matrices,vectors,trafo=[A1],[b2],[]
        go=spalte-1
        while go>-1:
            D1,d1,M1=build_M(A1,b2,go)
            matrices.append(D1)
            vectors.append(d1)
            trafo.append(M1)
            A1,b2=D1,d1
            go-=1
            
        #find first entry that is positive
            
        for i in range(len(d1)):
            if d1[i]>0:
                err=i
                break
        trafo.reverse()
        S=mul_list(trafo)[0]

        #multiply all the matrices: from left all the way to the right
        
        E=einheitsvektor(len(S),err)
        E=transpose(E)

        #create unit vector transposed with 1 on the err-th entry
        #multiply the transposed vector with S
        print((False,multiply(E,S)[0]))
        
        
        
                
        
            
            
            
            
            
        
        
        
        


if sys.argv[1]=="project":
    project(sys.argv[2],int(sys.argv[3]),sys.argv[4])

if sys.argv[1]=="image":
    image(sys.argv[2],sys.argv[3],sys.argv[4])

if sys.argv[1]=="H_representation":
    H_representation(sys.argv[2],sys.argv[3])

if sys.argv[1]=="compute_x_or_y":
    compute_x_or_y(sys.argv[2])




'''
For testing purposes

A1=[[1,1,2],[2,0,3],[1,0,-4],[-2,1,-1],[1,1,0]]
b1=[[2],[3],[4],[5],[1]]

A2=[[0.75, 0.5], [-1.5, 1.5], [0.9166666666666666, 0.0], [-1.3333333333333335, 1.0], [1, 1]]
b2=[[2.0], [6.0], [2.0], [6.0], [1]]
    

GM=[[0.0, -2.0, 2.0], [0.0, -1.0, 2.0], [1.0, -1.0, 1.0], [0.0, 2.0, -2.0], [0.0, 1.0, -2.0], [1.0, 1.0, -3.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, -1.0], [-1.0, -1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, 0.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, -2.0]]

    
GV=[[-1.0], [-2.0], [-2.0], [-2.0], [-1.0], [-2.0], [-2.0], [-2.0], [-1.0], [-1.0], [-1.0], [-1.0], [0.0], [-2.0], [-2.0]]
   
'''
    
        

    

    
    
    
            



    
    
    
                

            
            
        
            
    
    


        
        

        
