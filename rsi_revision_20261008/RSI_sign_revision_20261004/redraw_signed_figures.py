"""Publication figures from retained coordinates and initializer summaries."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib import font_manager
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "figures"
font_manager.findfont("Times New Roman", fallback_to_default=False)
plt.rcParams.update({"font.family": "Times New Roman", "font.size": 9.5,
    "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic", "mathtext.bf": "Times New Roman:bold",
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": .65})
BLUE, RED = "#245d9f", "#c24430"

def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def save(fig, name):
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight", pad_inches=.05)
    fig.savefig(OUT / (name + ".png"), bbox_inches="tight", pad_inches=.05, dpi=300)
    plt.close(fig)
    print("Saved", name)

def orbit():
    fig, ax = plt.subplots(figsize=(3.4,3.25))
    ax.set_aspect("equal"); ax.set(xlim=(-1.55,1.55),ylim=(-.52,2.55))
    beta=2*np.pi-.65; t=np.linspace(-np.pi/2,-np.pi/2+beta,400)
    xx, yy=np.cos(t),1+np.sin(t)
    ax.add_patch(Circle((0,1),1,fill=False,color=".7",ls="--",lw=.8))
    ax.plot(xx,yy,color=BLUE,lw=1.6)
    ax.annotate("",xy=(xx[190],yy[190]),xytext=(xx[175],yy[175]),
        arrowprops={"arrowstyle":"->","color":BLUE,"lw":1.5})
    ax.scatter([0],[0],s=28,color="#2f854f",zorder=4)
    ax.scatter([xx[-1]],[yy[-1]],s=30,color=RED,zorder=4)
    for endpoint in [(0,0),(xx[-1],yy[-1])]:
        ax.plot([0,endpoint[0]],[1,endpoint[1]],color=".45",ls="--",lw=.8,zorder=2)
    ax.scatter([0],[1],s=8,color=".35",zorder=4)
    ax.plot([0,xx[-1]],[0,yy[-1]],ls="--",color=RED,lw=1)
    ax.text(-.45,-.08,r"$\Delta$",color=RED)
    ax.annotate("Image-plane\nposition",xy=(xx[-1],yy[-1]),xytext=(-1.47,-.10),
        va="top",fontsize=9,arrowprops={"arrowstyle":"-","color":RED,"lw":.6})
    ax.text(.08,-.35,"Emission",color="#2f854f",fontsize=9)
    ax.plot([0,1],[1,1],color=".45",lw=.75)
    ax.text(.43,1.08,r"$r_L$")
    ax.plot(.38*np.cos(t),1+.38*np.sin(t),color="#237e80",lw=1)
    ax.add_patch(FancyArrowPatch((.38*np.cos(t[-12]),1+.38*np.sin(t[-12])),
        (.38*np.cos(t[-1]),1+.38*np.sin(t[-1])),arrowstyle="->",
        mutation_scale=9,color="#237e80",lw=1))
    ax.text(-.08,1.17,r"$\beta$",color="#237e80")
    ax.text(.37,2.32,r"$\mathbf{B}\;\odot$  (+z)",fontsize=10)
    ax.annotate("",xy=(.62,0),xytext=(.09,0),
        arrowprops={"arrowstyle":"->","color":BLUE,"lw":1})
    ax.text(.62,-.18,r"$\mathbf{v}_\perp$",color=BLUE)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); save(fig,"fig2_1_signed")

def spots():
    coords=rows(HERE/"case01_coordinates.csv")
    ids=np.array([int(r["pidx"]) for r in coords]); assert np.array_equal(ids,np.arange(1,50001))
    nr_tmm=np.array([[float(r["tmm_x_nm"]),float(r["tmm_y_nm"])] for r in coords])
    nr=np.array([[float(r["nr_full_x_nm"]),float(r["nr_full_y_nm"])] for r in coords])
    initial=ROOT/"RSI_revision_simulation_package_20260901"/"07_source_inputs"/"intrinsic_pair"/"01_intrinsic_error_initial.csv"
    a=np.genfromtxt(initial,delimiter=",",comments="%",encoding="utf-8-sig"); a=a[np.argsort(a[:,0])]
    assert np.array_equal(a[:,0],ids); launch=a[:,1:3]*1e3
    relfile=ROOT/"RSI大修仿真统一归档_20260911"/"03_COMSOL完整结果"/"revision_comsol"/"baseline_three_model_20260904"/"REL_full"/"poincare_zrel_1000um_with_pidx.csv"
    rel=np.genfromtxt(relfile,delimiter=",",comments="%",encoding="utf-8-sig"); rel=rel[np.argsort(rel[:,3])]
    assert np.array_equal(np.rint(rel[:,3]).astype(int),ids); rel=rel[:,:2]*1e3
    pick=np.sort(np.random.default_rng(2026100409).choice(50000,800,replace=False))
    fig,axes=plt.subplots(1,4,figsize=(7.1,2.6),sharex=True,sharey=True)
    circles=[("FW50","#2f854f","-"),("FW90","#b78112","--"),("FW99",BLUE,":")]
    for k,(ax,xy,name) in enumerate(zip(axes,[launch,nr_tmm,nr,rel],["Source","NR-TMM","NR-full","REL-full"])):
        center=xy.mean(axis=0); r=np.linalg.norm(xy-center,axis=1)
        ax.scatter(xy[pick,0],xy[pick,1],s=2,color=".3",alpha=.65,edgecolors="none",rasterized=False)
        if k>0:
            for q,(_,color,style) in zip([.5,.9,.99],circles):
                ax.add_patch(Circle(center,np.quantile(r,q,method="linear"),fill=False,color=color,ls=style,lw=1))
        else: ax.add_patch(Circle((0,0),18,fill=False,color=".5",lw=.7))
        ax.set_aspect("equal"); ax.set(xlim=(-23,23),ylim=(-23,23),xticks=[-20,0,20],yticks=[-20,0,20])
        ax.set_title(f"({chr(97+k)}) {name}"); ax.set_xlabel("x (nm)")
        ax.grid(color=".9",lw=.4)
    axes[0].set_ylabel("y (nm)")
    keys=[Line2D([],[],color=c,ls=s,lw=1,label=n) for n,c,s in circles]
    fig.legend(handles=keys,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.6,0))
    fig.subplots_adjust(left=.07,right=.995,top=.88,bottom=.25,wspace=.17)
    save(fig,"fig3_1_signed")

def sampling():
    bins=rows(ROOT/"RSI_reviewer8_20261004"/"sampling_histogram_bins.csv")
    gof=rows(ROOT/"RSI_reviewer8_20261004"/"sampling_gof.csv")
    fig,axes=plt.subplots(2,2,figsize=(7.1,4.8))
    for row,n in enumerate([50000,50000000]):
        for j,var in enumerate(["U0","alpha"]):
            ax=axes[row,j]; b=[v for v in bins if int(v["N"])==n and v["variable"]==var]
            g=next(v for v in gof if int(v["N"])==n and v["variable"]==var)
            low=np.array([float(v["bin_low"]) for v in b]); high=np.array([float(v["bin_high"]) for v in b])
            y=np.array([float(v["empirical_pdf"]) for v in b]); ci=np.array([[float(v["ci95_pdf_low"]),float(v["ci95_pdf_high"])] for v in b])
            bars=ax.bar((low+high)/2,y,width=high-low,color=BLUE,alpha=.4,label="Generated PDF")
            error=ax.errorbar((low+high)/2,y,yerr=np.maximum(np.vstack((y-ci[:,0],ci[:,1]-y)),0),fmt=".",color=".25",ms=1.5,lw=.55,capsize=1.5,label="Pointwise 95% CI")
            x=np.linspace(low[0],high[-1],400)
            theory=np.full_like(x,2) if j==0 else np.sin(2*np.deg2rad(x))/np.sin(np.deg2rad(10.7))**2*np.pi/180
            curve,=ax.plot(x,theory,color=RED,lw=1.3,label="Theoretical PDF")
            ax.set(xlim=(low[0],high[-1]),ylim=(0,2.5 if j==0 else .23))
            ax.set_xlabel(r"$K_0$ (eV)" if j==0 else r"Emission angle $\alpha$ (deg)")
            ax.set_ylabel("PDF (eV$^{-1}$)" if j==0 else "PDF (deg$^{-1}$)")
            d=r"$D_{\rm KS}$" if row==0 else r"$D_{\rm grid}$"
            title=f"({chr(97+row*2+j)}) " + ("Energy" if j==0 else "Angle")
            ns=r"$5\times10^4$" if row==0 else r"$5\times10^7$"
            ax.set_title(title+"; n = "+ns+"\n"+d+f" = {float(g['CDF_max_deviation']):.6f}")
            ax.grid(color=".9",lw=.4)
    fig.legend([bars,curve,error],["Generated PDF","Theoretical PDF","Pointwise 95% CI"],loc="upper center",ncol=3,frameon=False)
    fig.tight_layout(rect=(0,0,1,.93)); save(fig,"fig_sampling_reviewer8")

def stripes():
    data=rows(HERE/"stripe_paired_coordinates.csv"); metrics=rows(HERE/"stripe_metrics.csv")
    pick=np.array([int(v["pidx"])-1 for v in rows(HERE/"stripe_display_ids.csv")])
    assert len(pick)==3000 and len(np.unique(pick))==3000
    xvals=[float(v[key]) for v in data for key in ["tmm_x_nm","nr_full_x_nm"]]
    yvals=[float(v[key]) for v in data for key in ["tmm_y_nm","nr_full_y_nm"]]
    fig,axes=plt.subplots(2,1,figsize=(7.1,6.25))
    for k,(ax,name) in enumerate(zip(axes,["11_center","13_axis_minus_1um"])):
        c=[v for v in data if v["case_name"]==name]; assert [int(v["pidx"]) for v in c]==list(range(1,50001))
        xy=np.array([[float(v[key]) for key in ["tmm_x_nm","tmm_y_nm","nr_full_x_nm","nr_full_y_nm"]] for v in c])
        full=ax.plot(xy[pick,2],xy[pick,3],".",color=RED,ms=2,label="NR-full (COMSOL)")[0]
        tmm=ax.plot(xy[pick,0],xy[pick,1],"o",color=BLUE,ms=1.8,mfc="none",mew=.35,label="NR-TMM")[0]
        m={v["metric"]:float(v["difference"]) for v in metrics if v["case_name"]==name}
        name2="Centered magnetic axis" if k==0 else r"Magnetic axis x = $-1\,\mu$m"
        ax.set_title(f"({chr(97+k)}) "+name2+f"\n$d_c$ = {m['centroid_vector_error_nm']:.3f} nm; paired RMS = {m['paired_landing_RMS_nm']:.3f} nm")
        ax.set_aspect("equal"); ax.set(xlim=(min(xvals)-4,max(xvals)+4),ylim=(min(yvals)-4,max(yvals)+4),xlabel="x (nm)",ylabel="y (nm)")
        ax.grid(color=".9",lw=.4)
    fig.legend([tmm,full],["NR-TMM","NR-full (COMSOL)"],loc="upper center",ncol=2,frameon=False)
    fig.tight_layout(rect=(0,0,1,.95)); save(fig,"fig_stripes_reviewer8")

if __name__=="__main__":
    orbit(); spots(); sampling(); stripes()
