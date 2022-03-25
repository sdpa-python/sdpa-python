#!/usr/bin/env python
"""param.py

A module of sdpap
Parameter definition of sdpap

December 2010, Kenta KATO
"""

__all__ = ['param']

import multiprocessing # for cpu count

def param(option=None):
    """Create SDPAP parameters

    Args:
      option: A dictionary for parameters.
        example:
          {'epsilonStar':1.0e-7, 'epsilonDash':1.0e-5}

    Returns:
      option: A dictionary of parameters. Each field is as follows:
        * maxIteration: The maximum number of iterations.
        * epsilonStar: The accuracy of an approximate optimal solution
            for primal and dual SDP.
        * lambdaStar: An initial point.
        * omegaStar: The search region for an optimal solution.
        * lowerBound: Lower bound of the minimum objective value of
            the primal SDP.
        * upperBound: Upper bound of the maximum objective value of
            the dual SDP
        * betaStar: The parameter for controlling the search direction
            if the current point is feasible.
        * betaBar: The parameter for controlling the search direction
            if the current point is infeasible.
        * gammaStar: A reduction factor for the primal and dual step lengths.
        * epsilonDash: The relative accuracy of an approximate optimal
            solution between primal and dual SDP.
        * isSymmetric: The flag for the checking the symmetricity of input
            matrices. (False => no check, True => check)
        * isDimacs: The flag to compute DIMACS ERROR
            (False => no computation, True => computation)
        * xPrint: default %+8.3e, NOPRINT skips printout
        * yPrint: default %+8.3e, NOPRINT skips printout
        * sPrint: default %+8.3e, NOPRINT skips printout
        * infPrint: default %+10.16e, NOPRINT skips printout
        * print: Destination of file output.
            the default setting is stdout by 'display'.
            If print is set 'no' or empty, no message is print out
        * resultFile: Destination of detail file output
        * sdpaResult: Destination of file output for SDPA result
        * numThreads: Number of Threads for internal computation
        * convMethod: Conversion method from CLP to LMI or EQ standard form

        The following parameters are for Exploiting sparsity algorithm.
        * domainMethod: Algorithm option for domain space. Parameter is:
            'none': Exploiting no sparsity in the domain space
            'clique': Applying dconv_cliquetree
            'basis': Applying dconv_basisrep
        * rangeMethod: Algorithm option for range space. Parameter is:
            'none': Exploiting no sparsity in the range space
            'clique': Applying rconv_cliquetree
            'decomp': Applying rconv_matdecomp

        The following parameters are for Free variables elimination.
        * frvMethod: The method to eliminate free variables
                       ('split' or 'elimination')
        * rho: The parameter of range in split method or
                 pivoting in elimination method
        * zeroPoint: The zero point of matrix operation,
            determine unboundness, or LU decomposition


        One recommended choice of Exploiting sparsity parameters:
          domain = 'clique'
          range  = 'decomp'
          convMethod = 'EQ'
        The other recommended choice of parameters:
          domain = 'basis'
          range  = 'clique'
          convMethod = 'LMI'
    """

    # Default parameter
    option0 = {
        'maxIteration': 100,
        'epsilonStar': 1.0E-7,
        'lambdaStar': 1.0E2,
        'omegaStar': 2.0,
        'lowerBound': -1.0E5,
        'upperBound': 1.0E5,
        'betaStar': 0.1,
        'betaBar': 0.2,
        'gammaStar': 0.9,
        'epsilonDash': 1.0E-7,
        'isSymmetric': False,
        'isDimacs': False,
        'xPrint': '%+8.3e',
        'yPrint': '%+8.3e',
        'sPrint': '%+8.3e',
        'infPrint': '%+10.16e',
        'print': 'display',
        'resultFile': '',
        'sdpaResult': '',
        'numThreads': multiprocessing.cpu_count(), # Max Avialable Number
        'frvMethod': 'split',  # 'split' or 'elimination'
        'convMethod': 'LMI',    # 'LMI' or 'EQ'
        'domainMethod': 'none',  # 'none' or 'clique' or 'basis'
        'rangeMethod': 'none',  # 'none' or 'clique' or 'decomp'
        'rho': 0.0,
        'zeroPoint': 1.0E-12
        }

    if option == None:
        return option0
    if not isinstance(option, dict):
        print("Input parameter is not dictionary.\n"
              "Default parameter is set.")
        return option0

    if 'maxIteration' not in option:
        option['maxIteration'] = option0['maxIteration']
    elif not isinstance(option['maxIteration'], int):
        print("option.maxIteration must be integer.\n"
              "Default value is set.")
        option['maxIteration'] = option0['maxIteration']

    if 'epsilonStar' not in option:
        option['epsilonStar'] = option0['epsilonStar']
    elif not isinstance(option['epsilonStar'], float):
        print("option.epsilonStar must be float.\n"
              "Default value is set.")
        option['epsilonStar'] = option0['epsilonStar']

    if 'lambdaStar' not in option:
        option['lambdaStar'] = option0['lambdaStar'];
    elif not isinstance(option['lambdaStar'], float):
        print("option.lambdaStar must be float.\n"
              "Default value is set.")
        option['lambdaStar'] = option0['lambdaStar']

    if 'omegaStar' not in option:
        option['omegaStar'] = option0['omegaStar']
    elif not isinstance(option['omegaStar'], float):
        print("option.omegaStar must be float.\n"
              "Default value is set.")
        option['omegaStar'] = option0['omegaStar']

    if 'lowerBound' not in option:
        option['lowerBound'] = option0['lowerBound']
    elif not isinstance(option['lowerBound'], float):
        print("option.lowerBound must be float.\n"
              "Default value is set.")
        option['lowerBound'] = option0['lowerBound']

    if 'upperBound' not in option:
        option['upperBound'] = option0['upperBound']
    elif not isinstance(option['upperBound'], float):
        print("option. must be float.\n"
              "Default value is set.")
        option['upperBound'] = option0['upperBound']

    if 'betaStar' not in option:
        option['betaStar'] = option0['betaStar']
    elif not isinstance(option['betaStar'], float):
        print("option. must be float.\n"
              "Default value is set.")
        option['betaStar'] = option0['betaStar']

    if 'betaBar' not in option:
        option['betaBar'] = option0['betaBar']
    elif not isinstance(option['betaBar'], float):
        print("option.betaBar must be float.\n"
              "Default value is set.")
        option['betaBar'] = option0['betaBar']

    if 'gammaStar' not in option:
        option['gammaStar'] = option0['gammaStar']
    elif not isinstance(option['gammaStar'], float):
        print("option.gammaStar must be float.\n"
              "Default value is set.")
        option['gammaStar'] = option0['gammaStar']

    if 'epsilonDash' not in option:
        option['epsilonDash'] = option0['epsilonDash']
    elif not isinstance(option['epsilonDash'], float):
        print("option.epsilonDash must be float.\n"
              "Default value is set.")
        option['epsilonDash'] = option0['epsilonDash']

    if 'searchDir' in option:
        print("Parameter *searchDir* is no longer supported.\n"
              "HRVW/KSH/M is automatically used.")

    if 'isSymmetric' not in option:
        option['isSymmetric'] = option0['isSymmetric']
    elif not isinstance(option['isSymmetric'], bool):
        print("option.isSymmetric must be bool.\n"
              "Default value is set.")
        option['isSymmetric'] = option0['isSymmetric']

    if 'isDimacs' not in option:
        option['isDimacs'] = option0['isDimacs']
    elif not isinstance(option['isDimacs'], bool):
        print("option.isDimacs must be bool.\n"
              "Default value is set.")
        option['isDimacs'] = option0['isDimacs']

    if 'xPrint' not in option:
        option['xPrint'] = option0['xPrint']
    elif not isinstance(option['xPrint'], str):
        print("option.xPrint must be string.\n"
              "Default value is set.")
        option['xPrint'] = option0['xPrint']

    if 'yPrint' not in option:
        option['yPrint'] = option0['yPrint']
    elif not isinstance(option['yPrint'], str):
        print("option.YPrint must be string.\n"
              "Default value is set.")
        option['yPrint'] = option0['yPrint']

    if 'sPrint' not in option:
        option['sPrint'] = option0['sPrint']
    elif not isinstance(option['sPrint'], str):
        print("option.sPrint must be string.\n"
              "Default value is set.")
        option['sPrint'] = option0['sPrint']

    if 'infPrint' not in option:
        option['infPrint'] = option0['infPrint']
    elif not isinstance(option['infPrint'], str):
        print("option.infPrint must be string.\n"
              "Default value is set.")
        option['infPrint'] = option0['infPrint']

    if 'print' not in option:
        option['print']=option0['print']
    elif len(option['print']) == 0:
        option['print'] = 'no'
    elif not isinstance(option['print'], str):
        print("option.print must be string.\n"
              "Default value is set.\n"
              "*** option.print must be string for FILE. ***\n"
              "  \"display\" is for stdout.\n"
              "  \"no\" or empty is for skip message.\n"
              "  filename is filename in which message will be written.\n"
              "*** option.print must be string for FILE. ***")
        option['print'] = option0['print']

    if 'resultFile' not in option or len(option['resultFile']) == 0:
        option['resultFile'] = option0['resultFile']
    elif not isinstance(option['resultFile'], str):
        print("option.resultFile must be string.\n"
              "Default value is set.")
        option['resultFile'] = option0['resultFile']

    if 'sdpaResult' not in option or len(option['sdpaResult']) == 0:
        option['sdpaResult'] = option0['sdpaResult']
    elif not isinstance(option['sdpaResult'], str):
        print("option.sdpaResult must be string.\n"
              "Default value is set.")
        option['sdpaResult'] = option0['sdpaResult']

    if 'numThreads' not in option:
        option['numThreads'] = option0['numThreads']
    elif not isinstance(option['numThreads'], int):
        print("option.numThreads must be integer.\n"
              "Default value is set.")
        option['numThreads'] = option0['numThreads']

    if 'frvMethod' not in option:
        option['frvMethod'] = option0['frvMethod']
    elif len(option['frvMethod']) == 0:
        option['frvMethod'] = option0['frvMethod']
    elif not isinstance(option['frvMethod'], str):
        print("option.frvMethod must be 'split' or 'elimination'.\n"
              "Default value is set.")
        option['frvMethod'] = option0['frvMethod']
    elif (option['frvMethod'] != 'split' and
          option['frvMethod'] != 'elimination'):
        print("option.frvMethod must be \"split\" or \"elimination\".\n"
              "Default value is set.")
        option['frvMethod'] = option0['frvMethod']

    if 'convMethod' not in option:
        option['convMethod'] = option0['convMethod']
    elif len(option['convMethod']) == 0:
        option['convMethod'] = option0['convMethod']
    elif not isinstance(option['convMethod'], str):
        print("option.convMethod must be 'LMI' or 'EQ'.\n"
              "Default value is set.")
        option['convMethod'] = option0['convMethod']
    elif option['convMethod'] != 'LMI' and option['convMethod'] != 'EQ':
        print("option.toSedumi must be 'LMI' or 'EQ'.\n"
              "Default value is set.")
        option['convMethod'] = option0['convMethod']

    if 'domainMethod' not in option:
        option['domainMethod'] = option0['domainMethod']
    elif len(option['domainMethod']) == 0:
        option['domainMethod'] = option0['domainMethod']
    elif isinstance(option['domainMethod'], str) == False:
        print("option.domainMethod must be 'none' or 'clique' or 'basis'.\n"
              "Default parameter is set.")
        option['domainMethod'] = option0['domainMethod']
    elif option['domainMethod'] != 'none' and \
             option['domainMethod'] != 'clique' and \
             option['domainMethod'] != 'basis':
        print("option.domainMethod must be 'none' or 'clique' or 'basis'.\n"
              "Default parameter is set.")
        option['domainMethod'] = option0['domainMethod']

    if 'rangeMethod' not in option:
        option['rangeMethod'] = option0['rangeMethod']
    elif len(option['rangeMethod']) == 0:
        option['rangeMethod'] = option0['rangeMethod']
    elif isinstance(option['rangeMethod'], str) == False:
        print("option.rangeMethod must be 'none' or 'clique' or 'decomp'.\n"
              "Default parameter is set.")
        option['rangeMethod'] = option0['rangeMethod']
    elif option['rangeMethod'] != 'none' and \
             option['rangeMethod'] != 'clique' and \
             option['rangeMethod'] != 'decomp':
        print("option.rangeMethod must be 'none' or 'clique' or 'decomp'.\n"
              "Default parameter is set.")
        option['rangeMethod'] = option0['rangeMethod']

    if 'rho' not in option:
        option['rho'] = option0['rho']
    elif not isinstance(option['rho'], float):
        print("option.rho must be float.\n"
              "Default value is set.")
        option['rho'] = option0['rho']

    if 'zeroPoint' not in option:
        option['zeroPoint'] = option0['zeroPoint']
    elif not isinstance(option['zeroPoint'], float):
        print("option.zeroPoint must be float.\n"
              "Default value is set.")
        option['zeroPoint'] = option0['zeroPoint']

    return option
