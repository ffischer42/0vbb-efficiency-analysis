function get_share_for_worker(to_do_list, num_workers, this_worker)
    step = Int(round(length(to_do_list) / num_workers))
    if this_worker != num_workers
        set = ((this_worker - 1) * step +1):(this_worker*step)
    else
        set = ((this_worker - 1) * step +1):length(to_do_list)
    end
    return to_do_list[set]
end