#!/usr/bin/env bash
# Build htslib from source inside the cibuildwheel manylinux container.
#
# The _cppext extension links against libhts and libz. manylinux images ship
# neither htslib nor a usable htslib package, so it is compiled here and
# installed into $PREFIX; auditwheel then bundles libhts.so.3 into the wheel.
#
# Pinned to the same htslib the local verified build linked against, so the
# published wheels behave identically to a developer conda build.
set -euo pipefail

HTSLIB_VERSION="${HTSLIB_VERSION:-1.23.1}"
PREFIX="${HTSLIB_PREFIX:-/usr/local}"

echo "==> building htslib ${HTSLIB_VERSION} into ${PREFIX}"

# bzip2/xz are htslib's CRAM codecs. ScanMST only calls core BAM/SAM APIs, but
# they are cheap and keep htslib able to read CRAM if a user passes one.
if command -v dnf >/dev/null 2>&1; then
    dnf install -y zlib-devel bzip2-devel xz-devel
elif command -v yum >/dev/null 2>&1; then
    yum install -y zlib-devel bzip2-devel xz-devel
else
    echo "no dnf/yum available in this image" >&2
    exit 1
fi

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT
cd "${workdir}"

curl -fsSL -o htslib.tar.bz2 \
    "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
tar xf htslib.tar.bz2
cd "htslib-${HTSLIB_VERSION}"

# libcurl/gcs/s3 are deliberately off: ScanMST reads local files only, and
# enabling them would drag libcurl + openssl into the bundled wheel.
# libdeflate is off to avoid another bundled dependency.
./configure \
    --prefix="${PREFIX}" \
    --disable-libcurl \
    --disable-gcs \
    --disable-s3 \
    --disable-plugins \
    --without-libdeflate

make -j"$(nproc)"
make install

ldconfig || true

echo "==> htslib installed:"
ls -l "${PREFIX}/lib/libhts."* || true
test -f "${PREFIX}/include/htslib/hts.h" || { echo "missing headers" >&2; exit 1; }
