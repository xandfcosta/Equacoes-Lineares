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
        mat = list()
        i = 0
        for atual in linhas:
            atual = atual.split( " " )
            for j in range( len( atual ) ):
                atual[ j ] = float( atual[ j ] )
            
            atual.append( float( b[ i ] ) )
            i += 1
            mat.append( atual )
        
    return mat, n
        
def mostra_matriz( mat ):
    for linha in mat:
        print( "|\t", end="" )
        
        for num in linha:
            print( f"{num:.3f}\t", end="" )
            
        print( "|" )
        
    print()
          
def pivotiza( mat, pos_pivo ):
    maior = mat[ pos_pivo ][ pos_pivo ]
    pos_maior = pos_pivo
       
    for linha in range( len( mat ) ):
        if mat[ linha ][ pos_pivo ] > maior:
            maior = mat[ linha ][ pos_pivo ]
            pos_maior = linha
            
            
    mat[ pos_pivo ], mat[ pos_maior ] = mat[ pos_maior ], mat[ pos_pivo ]
            
    return mat, mat[ pos_pivo ][ pos_pivo ]
            
def gauss_jordan( mat, m ):  
    for i in range( m ):
        linha_pivo = mat[ i ]
        pivo = mat[ i ][ i ]
        if i == m and pivo == 0 or pivo == 0 and mat[ i ][ m ] == 0:
            break
        
        for j in range( m ):
            # if pivo == 0 and mat[ i ][ m ]:
            #     break
            
            if i != j:
                a = mat[ j ][ i ] / pivo
                for k in range( len( mat[ j ] ) ):
                    mat[ j ][ k ] -= a * linha_pivo[ k ]
                
    # Solução
    if pivo == 0 and mat[ m - 1 ][ m ] == 0:
        s = 'SPI'
    elif pivo == 0 or mat[ m - 1 ][ m ] == 0:
        s = 'SI'    
    else:
        s = list()
        for i in range( m ):
            s.append( mat[ i ][ m ] / mat[ i ][ i ] )
        
    
    return mat, s

def gauss_seidel( mat, m ):
    pass

def main():    
    for i in range( 4, 5 ):
        print( f"[!] Arquivo {i}" )
        
        mat, n = ler_arquivo( f"inputs/{i}.txt" )
        
        mat_res, s = gauss_jordan( mat, n )
        mostra_matriz(mat_res)
        print( f"S = {s}\n" )

          
if __name__ == '__main__':
    main()
        
        
        