
#ifdef WIN_R_BUILD
  #if __cplusplus < 201103L
    #define OLD_WIN_R_BUILD
  #else 
    #define NEW_WIN_R_BUILD
  #endif
#else
  #if __cplusplus < 201103L
    #error Error: ranger requires a real C++11 compiler. You probably have to update gcc. 
  #endif
#endif

