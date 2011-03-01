# -*- Autoconf -*-

## Autoconf macro for testing if compiler provides string with
## function names.
##

# CIT_HAVE_FUNTIONSTRING
#   Defines preprocessor macro __FUNCTION_NAME__.
# ------------
AC_DEFUN([CIT_FUNCTIONSTRING], [
  AC_LANG(C++)
  set_function_name=no

  AC_MSG_CHECKING([whether C++ compiler defines __PRETTY_FUNCTION__])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[]],
 	             [[const char* name = __PRETTY_FUNCTION__;]])],
    [AC_MSG_RESULT(yes)
     set_function_name=yes
     AC_DEFINE([__FUNCTION_NAME__], [__PRETTY_FUNCTION__], [Define __FUNCTION_NAME__ to __PRETTY_FUNCTION__.])],
    [AC_MSG_RESULT(no)])

  if test "$set_function_name" == no; then
    AC_MSG_CHECKING([whether C++ compiler defines __FUNCTION__])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[]],
 	               [[const char* name = __FUNCTION__;]])],
      [AC_MSG_RESULT(yes)
       set_function_name=yes
       AC_DEFINE([__FUNCTION_NAME__], [__FUNCTION__], [Define __FUNCTION_NAME__ to __FUNCTION__.])],
      [AC_MSG_RESULT(no)])
    fi

  if test "$set_function_name" == no; then
    AC_MSG_CHECKING([whether C++ compiler defines __func__])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[]],
 	               [[const char* name = __func__;]])],
      [AC_MSG_RESULT(yes)
       set_function_name=yes
       AC_DEFINE([__FUNCTION_NAME__], [__func__], [Define __FUNCTION_NAME__ to __func__.])],
      [AC_MSG_RESULT(no)])
    fi
]))


dnl end of file
