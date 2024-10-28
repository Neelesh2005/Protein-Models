import pandas as pd
import datetime
import numpy as np
import scipy.stats
def get_hfi_returns():
    hfi = pd.read_csv(
        "data/edhec-hedgefundindices.csv",
        header=0,
        index_col=0,
        parse_dates=True,
        na_values=-99.99,
    )
    hfi = hfi/100
    hfi.index = pd.to_datetime(hfi.index).to_period('M')
    return hfi

def get_ind_returns():
    ind = pd.read_csv("./data/ind30_m_vw_rets.csv",header=0,index_col=0,parse_dates=True)/100
    ind.index = pd.to_datetime(ind.index,format = "%Y%m").to_period('M')
    ind.columns = [column.strip() for column in ind.columns ]
    return ind
def drawdown(return_series:pd.Series):
    wealth_index = 1000*(1+return_series).cumprod()
    previous_peaks = wealth_index.cummax()
    drawdowns = (wealth_index-previous_peaks)/previous_peaks
    return pd.DataFrame({
            "Wealth":wealth_index,
            "Peaks":previous_peaks,
            "Drawdown":drawdowns
        })


def get_ffme_returns():
    date_parser = lambda x: datetime.strptime(x, '%Y%m')
    me_m = pd.read_csv(
        "data/Portfolios_Formed_on_ME_monthly_EW.csv",
        header=0,
        index_col=0,
        parse_dates=True,
        na_values=-99.99,
        date_format=date_parser
    )
    rets = me_m[['Lo 10','Hi 10']]
    rets.columns = ['Small Cap','Large Cap']
    rets = rets/100
    rets.index = pd.to_datetime(rets.index,format = "%Y%m").to_period('M')
    return rets
def semideviation(r):
    is_negative = r<0
    return r[is_negative].std(ddof=0)
def skewness(r):
    demeaned_r = r - r.mean()
    std = r.std(ddof=0)
    exp = (demeaned_r**3).mean()
    return exp/(std**3)


def kurtosis(r):
    demeaned_r = r - r.mean()
    std = r.std(ddof=0)
    exp = (demeaned_r**4).mean()
    return exp/(std**4)


def isNormal(r,level=0.01):
    statistic,p_value = scipy.stats.jarque_bera(r)
    return p_value>level
def annualized_ret(r,periods_per_year):
    comp_grow = (1+r).prod()
    n_periods = r.shape[0]
    return comp_grow**(periods_per_year/n_periods) -1 

def annulaized_vol(r,periods_per_year):
    return r.std(ddof=0)*(periods_per_year**0.5)

def sharpe_ratio(r,rf_rate,pds_per_year):
    rf_per_pd = (1+rf_rate)**(1/pds_per_year)-1
    excess_ret = r - rf_per_pd
    ann_ex_rt = annualized_ret(excess_ret,pds_per_year)
    ann_vol = annulaized_vol(r,pds_per_year)
    return ann_ex_rt/ann_vol
    
    
def portfolio_return(weights, returns):
    """
    Weights -> Returns
    """
    return weights.T @ returns

def portfolio_vol(weights, covmat):
    """
    Weights -> Vol
    """
    return (weights.T @ covmat @ weights)**0.5


def plot_eft(n_portf,er,cov,style='.-'):
    """
    Plots the 2-asset efficient frontier
    """
    weights = [np.array([w,1-w]) for w in np.linspace(0,1,n_portf)]
    rets = [portfolio_return(weights,er) for weights in weights]
    vols = [portfolio_vol(weights,cov) for weights in weights]
    ef = pd.DataFrame({
        "Returns":rets,
        "Volatility":vols
    })
    return ef.plot.line(x="Volatility",y="Returns",style=style)