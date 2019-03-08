    # if good_day is not None:
    #     print(good_day)

    #         if good_day is not None:
    #             break

    # day = good_day
    # init = good_init
    # step = main_step_3



    # good_day = None
    # good_init = None

    # for month in available_months:

    #     url = url_base + month
    #     req = requests.get(url)

    #     text = req.content.split('/</a>')

    #     available_days = []

    #     for t in text:
    #         if '>20' in t:
    #             available_days.append(t.split('>')[-1].split('>')[-1])

    #     for day in available_days:

    #         yyyy2, mm2, dd2 = int(day[0:4]), int(day[4:6]), int(day[6:8])

    #         url = url_base + month + '/' + day
    #         req = requests.get(url)

    #         text = req.content.split('</a>')

    #         available_inits = []

    #         for t in text:
    #             if 'grb2' in t and not 'align' in t:
    #             # if ('grb2' in t or 'inv' in t) and not 'align' in t:
    #                 # print(t)
    #                 available_inits.append(t.split('>')[-1].split('>')[-1])

    #         for init in available_inits:

    #             # print(init)

    #             yyyy2, mm2, dd2, hhhh2, hhh2 = int(init[6:10]), int(init[10:12]), int(init[12:14]), int(init[15:19])/100 , int(init[20:23])

    #             # print(request_time, init, yyyy2, mm2, dd2, hhhh2, hhh2)

    #             init_dt = dt.datetime(yyyy2, mm2, dd2, hhhh2)
    #             delta_t = request_time - init_dt

    #             # print(request_time, init_dt)

    #             if delta_t.seconds < 0 or delta_t.days < 0:
    #                 continue

    #             # 00 06 12 18 ... for ensemble members

    #             ens_step = '%02d' % (6*int(delta_t.seconds/(6.*60*60)))

    #             # 00 03 06 09 12 ... for 0.5 and 1.0 deg main runs

    #             main_step = '%02d' % (3*int(delta_t.seconds/(3.*60*60)))
    #             main_step_3 = '%03d' % (3*int(delta_t.seconds/(3.*60*60)))

    #             good_day = day
    #             good_init = init

    #             if good_day is not None:
    #                 break

    #         if good_day is not None:
    #             break

    # day = good_day
    # init = good_init
    # step = main_step_3