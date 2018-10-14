function opt = optsetfun(opt)

if ~exist('opt','var'),                              opt = [];      end;
if ~isfield(opt,'perplexity'),            opt.perplexity = 30;      end;
if ~isfield(opt,'epsilon'),                  opt.epsilon = 300;     end  % learning rate
if ~isfield(opt,'min_gain'),                opt.min_gain = 0.01;    end  % initial jitter
if ~isfield(opt,'momentum'),                opt.momentum = 0.5;     end  % initial momentum
if ~isfield(opt,'final_momentum'),    opt.final_momentum = 0.8;     end  % final momentum
if ~isfield(opt,'max_iter'),                opt.max_iter = 300;     end  % maximum number of iterations
if ~isfield(opt,'stop_lying_iter'),  opt.stop_lying_iter = 100;     end  % jitter decay
if ~isfield(opt,'mom_switch_iter'),  opt.mom_switch_iter = 250;     end  % iteration where momentum changes
if ~isfield(opt,'v'),                             opt.v = 1;        end  % iteration where momentum changes
if ~isfield(opt,'ReducedDim'),           opt.ReducedDim = 40;       end  % iteration where momentum changes