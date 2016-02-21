#include "quit_import.h"


bool is_open ( char c )
    {
    if ( c == '(' || c == '[' || c == '{' )
        {
        return true;
        }
    else
        {
        return false;
        }
    }

bool is_closed ( char c )
    {
    if ( c == ')' || c == ']' || c == '}' )
        {
        return true;
        }
    else
        {
        return false;
        }
    }

bool is_matching ( char cc, char co )
    {
    if ( ( cc == ')' && co == '(' ) or ( cc == ']' && co == '[' ) or ( cc == '}' && co == '{' ) )
        {
        return true;
        }
    else
        {
        return false;
        }
    }

void cleanup ( std::string  &text )
    {

    std::remove ( text.begin(),text.end(),'\n' );
    std::remove ( text.begin(),text.end(),'\r' );
    std::replace ( text.begin(),text.end(),'I',' ' );
    std::regex exp ( R"(\*10\^)" );
    text = std::regex_replace ( text,exp,"   e" );
    text.erase ( std::remove ( text.begin(),text.end(),' ' ),text.end() );

    }


std::complex< double > * parse_complex_matrix ( std::string text,int &size )
    {


    static std::complex< double > * matrix;

    std::regex cn ( R"([^(\r\,\{\}\s\n)][[:digit:]\^\+\-\.\*\s\n\rEexp]*I)" ); //recognieses mathematica complex numbers.

    std::regex rp ( R"([[:digit:]\+\-][[:digit:]\^\+\-\*\.\n\rEexp]*)" );

    std::regex ip ( R"([\+\-][\s\r\n]*[[:digit:]\^\+\-\*\.\\n\rEexp]*\sI$)" );

    auto M_begin = std::sregex_iterator ( text.begin(), text.end(), cn );
    auto M_end   = std::sregex_iterator();
    int s = sqrt ( std::distance ( M_begin, M_end ) );
    size = s;


    std::smatch d_match;

    std::cout
            << "Importierte Matrix hat Dimension:"
            << size
            << std::endl;

    matrix = new std::complex< double >[std::distance ( M_begin, M_end )];

    std::string d_match_str;

    int j = 0;
    for ( std::sregex_iterator i = M_begin ; i != M_end; ++i )
        {

        std::smatch match = *i;
        std::string match_str = match.str();

        std::regex_search ( match_str,d_match,rp );
        d_match_str = d_match.str();
        cleanup ( d_match_str );
        matrix[j].real ( atof ( d_match_str.c_str() ) );

        std::regex_search ( match_str,d_match,ip );
        d_match_str = d_match.str();
        cleanup ( d_match_str );
        matrix[j].imag ( atof ( d_match_str.c_str() ) );

        ++j;

        }
    return matrix;
    }

std::complex<double> * Mathematica_import_matrix_from_file ( std::string FILENAME, int &SIZE_IMPORTED )
    {

    std::ifstream FILE;
    FILE.open ( FILENAME );

    std::stack<char> parenthese;
    std::string text;
    std::string line;
    std::complex< double > * IMPORTED_MATRIX;

    if ( FILE.is_open() )
        {

        int * mat_count;

        while ( std::getline ( FILE,line ) )
            {
            text += line;
            text.push_back ( '\n' );
            }

        for ( int i = 0; i < text.length()-1; ++i )
            {
            if ( is_open ( text.at ( i ) ) )
                {
                parenthese.push ( text.at ( i ) );
                }
            if ( is_closed ( text.at ( i ) ) )
                {
                if ( parenthese.empty() )
                    {
                    std::cout << "WARNING!!! NOT MATCHING PARENTHESISES. Inconsitant Data" << std::endl;
                    }
                else
                    {
                    if ( is_matching ( text.at ( i ),parenthese.top() ) )
                        {
                            {
                            //std::cout << parenthese.top();
                            parenthese.pop();
                            }
                        }
                    else
                        {
                        std::cout << "WARNING!!! NOT MATCHING PARENTHESISES WILL TRY TO PARSE" << text.at ( i ) << parenthese.top() << std::endl;
                        }
                    }
                }


            }
        FILE.close();

        if ( parenthese.empty() )
            {
            IMPORTED_MATRIX = parse_complex_matrix ( text,SIZE_IMPORTED );
            return IMPORTED_MATRIX;
            }
        else
            {
            std::cout << "WILL TRY TO PARSE ANYWAY";
            IMPORTED_MATRIX = parse_complex_matrix ( text,SIZE_IMPORTED );
            return IMPORTED_MATRIX;

            }
        }
    else
        {

        std::cout << "Error Open File" << std::endl;

        }
    }
