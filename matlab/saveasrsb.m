function saveasrsb(fig,fn)
global outdir;

try
saveas(fig,[outdir '\' fn]);
catch
end