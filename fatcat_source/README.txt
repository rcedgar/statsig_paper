# To re-build patched FATCAT which reports
#	normalized score, p-value, mu, beta

git clone https://github.com/GodzikLab/FATCAT-dist.git
git checkout ef787fe

# Replace SigEva.C in FATCATMain/ with this patch
# Run make
