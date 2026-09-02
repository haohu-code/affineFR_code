function Robust_SAT_SDP()
% this is a wrapper for running SAT_SDP

load("robust_sdp_run.mat")

% initialize the OOM
output.objval = [];
output.runtime = [];
output.status = 'OOM';
save('robust_sdp_run.mat','output')
output = SAT_relaxation(model,relaxation);
save('robust_sdp_run.mat','output')

end