function [ tlocked, tstart ] = get_tlocked(pdz, t, events, params)

    tlocked = zeros(length(events),sum(params.prepost)+1);
    tstart = zeros(length(events),1);
    N=length(t);
        
    for i = 1 : length(events)

        idx = find(t>events(i),1,'first');
        if isempty(idx)
            warning('Event %d (%d) is beyond time vector.', i, events(i));
            break;
        end
        pdi = pdz(idx);
        pre = idx - params.prepost(1);
        if pre < 1
            pad0 = 1-pre;
            pre = 1;
        else
            pad0 = 0;
        end
        post = idx + params.prepost(2);
        if post > N
            pad1 = post-N;
            post = N;
        else
            pad1 = 0;
        end
        
        xx = pdz(pre:post);

        T = [nan(pad0,1);xx;nan(pad1,1)];
        tlocked(i,:) = T;
        tstart(i) = t(pre+pad0);

    end

end
