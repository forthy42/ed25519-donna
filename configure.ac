## ExampleLib Example: an example of using Automake to link with a library

AC_INIT([ED25519], [3:2], [bernd@net2o.de], [ed25519-donna],
        [http://keccak.noekeon.org/])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

CFLAGS=${CFLAGS--O3}

AC_PROG_CC
AM_PROG_AS
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

PRECC=${CC%gcc*}
PRECLANG=${RRECC%clang*}
LIBTOOL=${PRECLANG}libtool
AC_SUBST(LIBTOOL)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT