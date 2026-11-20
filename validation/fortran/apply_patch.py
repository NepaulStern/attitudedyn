"""Write alecia_oracle.f and the method=0 parameter file used by the oracle."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("fetch_vendor", ROOT / "fetch_vendor.py")
if _spec is None or _spec.loader is None:
    raise SystemExit("Could not load fetch_vendor.py")
_fetch = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fetch)
_fetch.ensure_downloaded(ROOT / "vendor")
SRC = (ROOT / "vendor" / "alecia_parallel.f").read_text()

MAIN_OLD = """      if(method.eq.0)then
         call parset0(parm,nparm)
         if(model.eq.1)call likely11(parm,nparm,obj,*3)
         if(model.eq.2)call likely12(parm,nparm,obj,*3)
   3     continue
         stop
         endif
"""
MAIN_NEW = """      if(method.eq.0)then
         call parset0(parm,nparm)
         if(model.eq.1)call likely11(parm,nparm,obj,*3)
         if(model.eq.2)call likely12(parm,nparm,obj,*3)
   3     continue
         call parset212(penalty)
         write(6,901)obj,penalty
 901     format(1x,'ORACLE_OBJ ',d25.16,1x,d25.16)
         open(22,file='oracle_obj.txt',status='replace')
         write(22,902)obj,penalty
 902     format(d25.16,1x,d25.16)
         close(22)
         call oracle_dump
         stop
         endif
"""

ANORM_OLD = """      subroutine anormc(arg,val)
c 
c       this subroutine evaluates the standard normal distribution
c       function.
c 
      real xt,valt
      real*8 arg,val
c 
c       standardize the argument of the distribution function 
c 
      xt=sngl(arg)*.707107
"""
ANORM_NEW = """      subroutine anormc(arg,val)
c 
c       this subroutine evaluates the standard normal distribution
c       function. ATTDYN_SHARED_MATH=1 uses the package's double Phi.
c 
      real xt,valt
      real*8 arg,val
      integer ishared
      call oracle_shared(ishared)
      if(ishared.eq.1)then
         call attdyn_cdf(arg,val)
         return
         endif
c 
c       standardize the argument of the distribution function 
c 
      xt=sngl(arg)*.707107
"""

# Only normadj1, whose next executable line is a bare return.
NORM_OLD = """      utc=(fac1*fuk)+pfac1
      utcr=sngl(utc)
      call utrunc(utcr)
      call mdnris(utcr,etcr,ier)
      utc=dble(etcr)
      return
      end
"""
NORM_NEW = """      utc=(fac1*fuk)+pfac1
      call oracle_shared(ishared)
      if(ishared.eq.1)then
         if(utc.gt.0.9999d0)utc=0.9999d0
         if(utc.lt.0.0001d0)utc=0.0001d0
         call attdyn_qnorm(utc,utc)
         return
         endif
      utcr=sngl(utc)
      call utrunc(utcr)
      call mdnris(utcr,etcr,ier)
      utc=dble(etcr)
      return
      end
"""

LIKE_OLD = """      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
"""
# This common header appears in several routines. Patch only likely12 by
# including the two lines that follow it there.
LIKE_INSERT_AFTER = """      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      call parset1(parm)
      call parset212(penalty)
"""
LIKE_INSERT = """      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/oracleft/ftall(nobss)
      call parset1(parm)
      call parset212(penalty)
"""

LOOP_OLD = """      f=0.d0
!$omp parallel do private (i,ft) reduction (+:f)
      do 2 i=1,nobs
         if(incluv(i).eq.1)then
            call obsmake2(i)
            call likobs12(ft)
            f=f+ft
            endif
   2     continue
!$omp end parallel do
"""
LOOP_NEW = """      f=0.d0
      do 2 i=1,nobs
         ft=0.d0
         if(incluv(i).eq.1)then
            call obsmake2(i)
            call likobs12(ft)
            f=f+ft
            endif
         ftall(i)=ft
   2     continue
"""


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{label}: expected 1 match, found {count}")
    return text.replace(old, new, 1)


def patch_likely12(text):
    start = text.find("      subroutine likely12(parm,nparm,f,*)")
    end = text.find("      subroutine likely2(", start)
    if start < 0 or end < 0:
        raise SystemExit("could not isolate likely12")
    body = text[start:end]
    if body.count(LIKE_INSERT_AFTER) != 1:
        raise SystemExit("likely12 progparms common was not unique")
    if body.count(LOOP_OLD) != 1:
        raise SystemExit("likely12 loop was not unique")
    body = body.replace(LIKE_INSERT_AFTER, LIKE_INSERT, 1)
    body = body.replace(LOOP_OLD, LOOP_NEW, 1)
    return text[:start] + body + text[end:]


text = SRC
text = replace_once(text, MAIN_OLD, MAIN_NEW, "main")
text = replace_once(text, ANORM_OLD, ANORM_NEW, "anormc")
text = replace_once(text, NORM_OLD, NORM_NEW, "normadj1")
text = patch_likely12(text)
(ROOT / "build").mkdir(exist_ok=True)
(ROOT / "build" / "alecia_oracle.f").write_text(text)
print("wrote", ROOT / "build" / "alecia_oracle.f")


def fortran_values(path):
    """List-directed values from alecia.p, stopping at the first label."""
    tokens = []
    lines = path.read_text().splitlines()
    label_at = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped[:1].isalpha():
            label_at = i
            break
        tokens.extend(stripped.split())
    header = tokens[:16]
    nseg = [int(float(x)) for x in tokens[16:22]]
    rest = tokens[22:]
    ncoef = (
        6 * 8
        + 6 * 10
        + 6 * (2 + 6)
        + (6 * 5 // 2)
        + sum(k - 2 for k in nseg)
    )
    # ncoef formula in parget for model 2:
    # (ny*ntime)+(ny*nbeta)+(ny*(nfac+ny))+(ny*(ny-1)/2)+sum(nseg-2)
    ny, nx, ntime, nfac = 6, 10, 8, 2
    ncoef = (ny * ntime) + (ny * nx) + (ny * (nfac + ny)) + (ny * (ny - 1) // 2)
    ncoef += sum(k - 2 for k in nseg)
    coef = [float(x) for x in rest[:ncoef]]
    icoef = [int(float(x)) for x in rest[ncoef : ncoef + ncoef]]
    if len(coef) != ncoef or len(icoef) != ncoef:
        raise SystemExit(f"parsed {len(coef)} coefs and {len(icoef)} flags, expected {ncoef}")
    labels = []
    for line in lines[label_at:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("method"):
            break
        labels.append(stripped)
    if len(labels) < ncoef:
        raise SystemExit(f"only {len(labels)} labels")
    labels = labels[:ncoef]
    return header, nseg, coef, icoef, labels


def x_vector(log_path):
    lines = log_path.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if "iteration = 8" in line and "obj" in line:
            start = i
            break
    if start is None:
        raise SystemExit("iteration 8 not found")
    block = []
    for line in reversed(lines[:start]):
        if "X VECTOR" in line:
            break
        block.append(line)
    block.reverse()
    names = []
    values = []
    for line in block:
        parts = line.split()
        if len(parts) != 2:
            continue
        name, raw = parts
        values.append(float(raw.replace("D", "E")))
        names.append(name)
    return names, values


header, nseg, coef, icoef, labels = fortran_values(ROOT / "src" / "alecia.p.gss")
names, values = x_vector(ROOT / "vendor" / "gss_log.txt")
free_labels = [labels[i] for i in range(len(labels)) if icoef[i] == 0]
if free_labels != names:
    for i, (a, b) in enumerate(zip(free_labels, names)):
        if a != b:
            raise SystemExit(f"label mismatch at free {i}: file {a} log {b}")
    raise SystemExit(f"free count file {len(free_labels)} log {len(names)}")
free_slot = 0
for i, flag in enumerate(icoef):
    if flag == 0:
        coef[i] = values[free_slot]
        free_slot += 1
header = list(header)
header[0] = "0"  # method = evaluate once
out = ROOT / "build" / "alecia.p"
chunks = []
chunks.append(" ".join(header))
chunks.append(" ".join(str(k) for k in nseg))


def emit(nums, per_line, fmt):
    lines = []
    for i in range(0, len(nums), per_line):
        lines.append(" ".join(fmt.format(x) for x in nums[i : i + per_line]))
    return lines


chunks.extend(emit(coef, 5, "{:.16g}"))
chunks.extend(emit(icoef, 5, "{:d}"))
chunks.append("")
chunks.extend(labels)
out.write_text("\n".join(chunks) + "\n")
print(f"spliced {free_slot} free coefficients into {out}")
