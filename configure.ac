## ExampleLib Example: an example of using Automake to link with a library

AC_INIT([ED25519], [3:13], [bernd@net2o.de], [ed25519-donna],
        [http://ed25519.cr.yp.to/])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AM_PROG_AS
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT

CFLAGS=`echo "${CFLAGS}" | sed -e 's/-O2//g'`

flags="$CFLAGS"
nomflags="`echo $CFLAGS | sed -e 's/-m[0-9]* //g'`"

# special optimizations for some plattforms
case "$host" in
    arm*-*linux*)
	flags="$CFLAGS"
	fastflags="-O3 -mfpu=neon $nomflags"
	enable_fast=yes
	CFLAGS="-O3 $nomflags"
	;;
    i686-*linux*|i686-*cygwin*)
	flags="$CFLAGS"
	fastflags="-O3 -m32 -mssse3 $nomflags"
	enable_fast=yes
	CFLAGS="-O2 -m32 $nomflags"
	;;
    i686-*)
	# if you don't want sse, specify i586-linux-gnu
	CFLAGS="-O2 -m32 -mssse3 $nomflags"
	;;
    *)
	# at least use O3, otherwise we might not get conditional moves
	CFLAGS="-O3 $nomflags"
	;;
esac

AC_ARG_ENABLE([fast],
[  --enable-fast    Turn on optional fast library with non-standard CPU features],
[case "${enableval}" in
  yes) fast=true ;;
  no)  fast=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-fast]) ;;
esac],[fast=false])
AM_CONDITIONAL([FAST], [test x$fast = xtrue])

PRECC=${CC%gcc*}
PRECLANG=${PRECC%clang*}
dnl MacOS X has a libtool that does something else
AC_CHECK_PROGS(LIBTOOL,"glibtool --tag=CC" libtool)
#LIBTOOL=${PRECLANG}${LIBTOOL}
echo "Libtool: $PRECC, $PRECLANG, $LIBTOOL"

AC_CHECK_FUNCS(explicit_bzero)
AC_SUBST(fastflags)
AC_SUBST(flags)
AC_SUBST(LIBTOOL)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
