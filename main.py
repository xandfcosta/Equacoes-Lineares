import copy
from decimal import *

MAX_FILES = 13

def ler_arquivo( caminho ):
    with open( caminho, "r" ) as arq:
        linhas = arq.readlines()
        
        # n = primeira linha
        n = int( linhas[ 0 ] )
        linhas.pop( 0 )
        
        # b[] = última linha
        b = linhas[ len( linhas ) - 1 ]
        linhas.pop( len( linhas ) - 1 )
        b = b.split( " " )
        
        # Organiza os dados do arquivo em uma matriz
        mat = []
        for i in range( len( linhas ) ):
            linhas[ i ] = linhas[ i ].split( " " )
            
            for j in range( len( linhas[ i ] ) ):
                linhas[ i ][ j ] = float( linhas[ i ][ j ] )
            linhas[ i ].append( float( b[ i ] ) )
            
            mat.append( linhas[ i ] )
        
    return mat, n
        
def mostra_matriz( mat ):
    for linha in mat:
        print( "|", end="" )
        
        for num in linha:
            if num == 0:
                num = abs( num ) # Arruma o -0
                
            print( f"{ format( num, '.3f' ):^10}", end="" )
            
        print( "|" )
        
    print()
        
# ============        
# GAUSS_JORDAN        
# ============
          
def pivotiza( mat, pos_pivo, m ):
    # Maior inicia como sendo o pivo
    maior = mat[ pos_pivo ][ pos_pivo ]
    pos_maior = pos_pivo
       
    # Verifica da linha do pivo pra frente
    for linha in range( pos_pivo, m ):
        if abs( mat[ linha ][ pos_pivo ] ) > maior:
            maior = abs( mat[ linha ][ pos_pivo ] )
            pos_maior = linha
            
    # Faz o swap das linhas
    mat[ pos_pivo ], mat[ pos_maior ] = mat[ pos_maior ], mat[ pos_pivo ]
            
def gauss_jordan( mat, m ):
    # Inicialização pra usar a biblioteca Decimal
    getcontext( )
    getcontext( ).prec = 50 # Precisão
    
    for j in range( m ): # Linha pivo
        pivotiza( mat, j, m )
        pivo = mat[ j ][ j ]
        
        # Caso ache um erro retorna imdediatamente
        if pivo == 0 and mat[ i ][ m ] == 0:
            return mat, "SPI"    
        elif pivo == 0 and mat[ i ][ m ] != 0:
            return mat, "SI"
        
        mat[ j ] = [ Decimal( x ) / Decimal( pivo ) for x in mat[ j ] ] # Divide a linha toda pelo pivo
        v = copy.deepcopy( mat[ j ] )
        
        for i in range( m ): # Outras linhas
            if i != j:
                a = Decimal( mat[ i ][ j ] ) / Decimal( mat[ j ][ j ] )
                
                for col in range( m + 1 ):
                    mat[ j ][ col ] = Decimal( mat[ j ][ col ] ) * Decimal( a )                 # Lj <- Lj * aij
                    mat[ i ][ col ] = Decimal( mat[ i ][ col ] ) - Decimal( mat[ j ][ col ] )   # Li <- Li - Lj
                    mat[ j ][ col ] = v[ col ]                                                  # Lj <- V
                
    # Solução   
    s = []
    for i in range( m ):
        s.append( float( mat[ i ][ m ] ) )
    
    
    return mat, s

# ============
# GAUSS-SEIDEL
# ============

def equaciona( mat, x, pos_linha, m ):
    soma = 0
        
    for i in range( m ):        
        if i != pos_linha:
            soma += mat[ pos_linha ][ i ] * x[ i ] # Faz a somatória dos outros x direto em uma variável
            
    return ( mat[ pos_linha ][ m ] - soma ) / mat[ pos_linha ][ pos_linha ] # b - ( tudo menos o x atual ) / x

def criterio_sassenfeld( mat, m ):
    # Vetor para empregar os Xi já obtidos
    # Inicia-se com 1 para não mudar nada na primeira iteração
    a = [ 1 for i in range( m ) ]
    
    for i in range( m ): # Linhas
        soma = 0
        pivotiza( mat, i, m )
        
        for j in range( m ): # Colunas
            if i != j:
                # Multiplica o absoluto do número pelo seu respectivo X anterior e adiciona na soma
                soma += abs( mat[ i ][ j ] ) * a[ j ]

        a[ i ] = ( soma / abs( mat[ i ][ i ] ) )
        
        # Caso seja maior ou igual a 1 já é inconclusivo
        if a[ i ] >= 1:
            return "Inconclusivo"
        
    return "Converge"
        
def gauss_seidel( mat, m, k, e ):
    x_list = [ 0 for i in range( m ) ]
    x_ant = copy.deepcopy( x_list )
    criterio = criterio_sassenfeld( mat, m )
    
    # Mostra as equações
    for i in range( m ):
        print( f"x{i + 1} = ( {mat[ i ][ m ]} - [ ", end="" ) 
        
        for j in range( m ):        
            if j != i:
                print( f"( {mat[ j ][ i ]}x{j + 1} ) {'+ ' if i < m - 1 else ''}",end="" )
                
        print( f") / {mat[ i ][ i ]}" )
        
        # Como divisão por 0 não existe logo o sistema é impossível
        if mat[ i ][ i ] == 0:
            return criterio, "SI"
    
    # Iterações
    for i in range( k ): # K
        print( f"k = {i + 1}" )
        cont = m
        
        for j in range( m ): # X   
            x = equaciona( mat, x_list, j, m )
            x_list[ j ] = x
            diferenca = abs( x - x_ant[ j ] )
            
            print( f"x{j + 1} = {x} | {diferenca} < e = {e} ", end="" ) 
            
            if diferenca < e:
                print( "( V )\n" )
                cont -= 1 # Controle para saber se todo X é válido
            else:
                print( "( F )\n" )
            
        print( "\n" )
        
        # Se todos X são válidos
        if cont == 0:
            # Retorna os valores com o critério
            return criterio, x_list
        else:
            # Senão copia os valores de X para o vetor de anterior e continua o processo
            x_ant = copy.deepcopy( x_list )

    # Caso tenha estourado o K retorna o que tiver
    return criterio, x_list

def main():    
    k = 1000
    e = 0.001
    
    for i in range( MAX_FILES ):
        print( f"[!] Arquivo {i}" )
        
        mat, n = ler_arquivo( f"inputs/{i}.txt" )
        
        mat_res, s = gauss_jordan( copy.deepcopy( mat ), n )
        
        print( "GAUSS-JORDAN" )
        mostra_matriz( mat_res )
        print( f"S = {s}\n" )
        
        print( "="*50 )
        
        print( "GAUSS-SEIDEL" )
        criterio, s = gauss_seidel( copy.copy( mat ), n, k, e )  
        print( f"[!] Criterio: {criterio}" )      
        print( f"S = {s}\n" )

          
if __name__ == '__main__':
    main()
        
        
        