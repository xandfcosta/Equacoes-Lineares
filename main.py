import copy
from os import system
from decimal import *

MAX_FILES = 13

def ler_arquivo( caminho ) -> tuple[ list, int ]:
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
                linhas[ i ][ j ] = Decimal( ( linhas[ i ][ j ] ) )
            linhas[ i ].append( Decimal( ( b[ i ] ) ) )
            
            mat.append( linhas[ i ] )
        
    return mat, n
        
def mostra_matriz( mat, pos = None ):
    i = 0
    
    for linha in mat:
        if i == pos:
            print( "\033[91m", end="" ) # Destaca a linha do erro com a cor vermelha
            
        print( "|", end="" )
        for num in linha:
            aux = str( num ).split( '.' )
            if aux[ 0 ] == '-0':
                num = 0.0
                
            print( f"{ format( num, '.3f' ):^13}", end="" )
            
        print( "|\033[00m" )            # Volta para a cor normal
        i += 1
        
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
        if mat[ linha ][ pos_pivo ].copy_abs( ) > maior:
            maior = mat[ linha ][ pos_pivo ].copy_abs( )
            pos_maior = linha
            
    # Faz o swap das linhas
    mat[ pos_pivo ], mat[ pos_maior ] = mat[ pos_maior ], mat[ pos_pivo ]
            
def gauss_jordan( mat, m ):
     
    # Linha pivo
    for i in range( m ): 
        pivotiza( mat, i, m )
        pivo = mat[ i ][ i ]
        
        # Caso ache um erro retorna imdediatamente
        # Retorna Matriz, solução e posição da linha atual
        if pivo == 0 and mat[ i ][ m ] == 0:
            return mat, "SPI", i                                    
        elif pivo == 0 and mat[ i ][ m ] != 0:
            return mat, "SI", i
        
        mat[ i ] = [ x / pivo for x in mat[ i ] ]  # Divide a linha toda pelo pivo
        v = copy.deepcopy( mat[ i ] )
        
        # Outras linhas
        for j in range( m ): 
            if i != j:
                a = mat[ j ][ i ] / mat[ i ][ i ]
                
                for col in range( m + 1 ):
                    mat[ i ][ col ] *= a                  # Li <- Li * aij
                    mat[ j ][ col ] -= mat[ i ][ col ]    # Lj <- Lj - Li
                    mat[ i ][ col ] = v[ col ]            # Li <- V
                
    # Solução   
    s = []
    for i in range( m ):
        s.append( float( mat[ i ][ m ] ) )
    
    return mat, s, None

# ============
# GAUSS-SEIDEL
# ============

def equaciona( mat, x, pos_linha, m ):
    soma = Decimal( '0' )
        
    # Faz a somatória dos outros x direto em uma variável
    for i in range( m ):        
        if i != pos_linha:
            soma += mat[ pos_linha ][ i ] * x[ i ]
            
    return ( mat[ pos_linha ][ m ] - soma ) / mat[ pos_linha ][ pos_linha ] # b - ( tudo menos o x atual ) / x

def troca_pivo( mat, pos, m ) -> None:
    pivo = mat[ pos ][ pos ]
    l_pivo = pos
    c_pivo = pos
    
    for l in range( pos, m ):
        for c in range( pos, m ):
            if mat[ l ][ c ].copy_abs( ) > pivo.copy_abs( ):
                pivo = mat[ l ][ c ]
                l_pivo = l
                c_pivo = c
            
    if l_pivo != pos:
        mat[ pos ], mat[ l_pivo ] = mat[ l_pivo ], mat[ pos ]
            
    if c_pivo != pos:
        for l in range( m ):
            mat[ l ][ pos ], mat[ l ][ c_pivo ] = mat[ l ][ c_pivo ], mat[ l ][ pos ]

def pivotamento_completo( mat, pos, m ) -> tuple[ str | None, int | None ]:
    for p in range( pos, m ):
        troca_pivo( mat, p, m )
        
        if mat[ p ][ p ].copy_abs( ) == Decimal( '0' ):
            if p != m - 1:
                troca_pivo( mat, p, m )
            else:
                return ( "SPI", p ) if mat[ p ][ m ].copy_abs( ) == Decimal( '0' ) else ( "SI", p )
        
        # Zera números abaixo do pivo
        v = copy.deepcopy( mat[ p ] )
        
        for l in range( p + 1, m ):
            a = mat[ l ][ p ] / mat[ p ][ p ]
            
            for c in range( m + 1 ):
                mat[ p ][ c ] *= a                # Li <- Li * aij
                mat[ l ][ c ] -= mat[ p ][ c ]    # Lj <- Lj - Li
                mat[ p ][ c ] = v[ c ]            # Li <- V
                
    return None, None

def criterio_sassenfeld( mat, m ) -> str:
    a = [ Decimal( '1' ) for _ in range( m ) ]
    maior = 0
    
    # Linhas
    for i in range( m ): 
        soma = Decimal( '0' )
        
        # Colunas
        for j in range( m ):
            if i != j:                          
                soma += mat[ i ][ j ].copy_abs( ) * a[ j ]

        a[ i ] = soma / mat[ i ][ i ].copy_abs( )
        if a[ i ] > maior:
            maior = a[ i ]
 
    return "Converge" if maior < 1 else "Inconclusivo"
                   
def gauss_seidel( mat, m, k, e ) -> tuple[ str | None, list | str, int | None ]:
    x_list = [ Decimal( '0' ) for _ in range( m ) ]
    x_ant = copy.deepcopy( x_list )
    criterio = None
    erro, linha_erro = pivotamento_completo( mat, 0, m )   
    
    if erro:
        return None, erro, linha_erro
    
    mostra_matriz(mat)
    criterio = criterio_sassenfeld( mat, m )
     
    # Mostra a equação do pivo atual
    for l in range( m ):
        print( f"x{l + 1} = ( {mat[ l ][ m ]:.10f} - [ ", end="" ) 
        for c in range( m ):                     
            if c != l:
                print( f"( {mat[ l ][ c ]:.10f} x {c + 1} ) {'+ ' if c < m - 1 else ''}", end="" )    
        print( f") / {mat[ l ][ l ]:.10f}" )
         
    # Iterações
    for i in range( k ): # K
        print( f"k = {i + 1}" )
        cont = m
        
        for j in range( m ): # X   
            x = equaciona( mat, x_list, j, m )
            x_list[ j ] = x
            diferenca = ( ( x ) - ( x_ant[ j ] ) ).copy_abs( )
            
            print( f"x{j + 1} = {x:.10f} | {diferenca:.10f} < e = {e} ", end="" ) 
            print( "( V )" if diferenca < e else "( F )" )
            
            if diferenca < e:
                cont -= 1 # Controle para saber se todo X é válido
            
        print( "\n" )
        
        # Se todos X são válidos
        if cont == 0:
            # Retorna os valores com o critério
            return criterio, x_list, linha_erro
        else:
            # Senão copia os valores de X para o vetor de anterior e continua o processo
            x_ant = copy.deepcopy( x_list )

    # Caso tenha estourado o K retorna o que tiver
    return criterio, x_list, linha_erro

def main(): 
    getcontext().prec = 40   
    k = 100
    e = 0.01
    
    for i in range( 9, 10 ):
        print( f"[!] Arquivo {i}" )
        
        mat, n = ler_arquivo( f"inputs/{i}.txt" )     
        
        print( "GAUSS-JORDAN" )
        mat_res, solucao, linha = gauss_jordan( copy.deepcopy( mat ), n )
        mostra_matriz( mat_res, linha )
        
        if linha != None:
            print( f"S = {solucao}" )
        else:      
            print( "S = { ", end= "" )
            for num in solucao:
                print( f"{num:.10f}, ", end= "")
            print( "}" )
        
        print( "-"* ( 10 * n + 2 ), "\n" )
        
        print( "GAUSS-SEIDEL" )
        criterio, solucao, linha = gauss_seidel( mat, n, k, e )  
        
        if linha != None:
            mostra_matriz( mat, linha )
            print( f"S = {solucao}" )
        else:
            print( f"[!] Criterio: {criterio}" )      
            
            print( "S = { ", end= "" )
            for num in solucao:
                print( f"{num:.10f}, ", end= "")
            print( "}" )

        print( "="* ( 10 * n + 2 ), "\n" )
        
        input()
        system("cls")
          
if __name__ == '__main__':
    main()
        
        
        