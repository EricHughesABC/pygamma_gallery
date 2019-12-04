from IPython.display import display, Math, Latex
import pygamma as pg

def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)

def bmatrix(aaa):
    """Returns a LaTeX bmatrix

    :aaa: numpy array or pg.matrix
    :returns: LaTeX bmatrix as a string
    """
    output = "\\begin{bmatrix} "

    nrows,ncols = aaa.shape
    for r in range(nrows):
        for c in range(ncols):
            realp = aaa[r][c].real
            imagp = aaa[r][c].imag
            
            if abs(realp) < 1e-16:
                realp=0.0
            if abs(imagp) < 1e-16:
                imagp=0.0

            if imagp == 0.0:
                str_num = "{:.4}".format(realp)
            elif realp == 0.0:
                str_num = "{:.4}j".format(imagp)
            else:
                str_num = "{:.4}".format(aaa[r][c])

            output += str_num


            if c+1 ==ncols:
                output +=  " \\\\ "
            else:
                output +=  " & "

    output = output + "\\end{bmatrix}"
    
    return output

def printLatex(a):
    
    if isinstance( a, pg.matrix):
        display(Latex(bmatrix(a.toNParray())))
    else:
        display(Latex(bmatrix(a)))
