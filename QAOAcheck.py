import qiskit
from qiskit import *
from qiskit.transpiler import CouplingMap
import networkx as nx
import numpy as np
import argparse
from fileinput import filename
from qiskit.visualization.pulse_v2 import draw
from qiskit.visualization import dag_drawer
from qiskit.converters import circuit_to_dag
from qiskit.transpiler import CouplingMap
from qiskit.visualization import plot_histogram
from qiskit.providers.aer.noise import NoiseModel
import re
import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import XGate
from qiskit.transpiler import PassManager, InstructionDurations
from qiskit.transpiler import PassManager, InstructionDurations

# from qiskit_ibm_provider.transpiler.passes.scheduling import DynamicCircuitScheduleAnalysis

# from qiskit.transpiler.passes import ALAPScheduleAnalysis
# from qiskit.visualization import timeline_drawer
from qiskit.transpiler import InstructionDurations
from os import walk
import os
import csv


def update_map(p2l, i, j):
    temp = p2l[i]
    p2l[i] = p2l[j]
    p2l[j] = temp
    return


def get_2D_layout(GG, row, column, cm):
    """_summary_

    Args:
        GG (_type_): Graph
        num (_type_): Number of nodes

    Returns:
        _type_: layout information using sabre
    """
    qc = QuantumCircuit(row * column)

    for e in GG.edges():
        qc.cx(e[0], e[1])
    cpp = cm
    temp = []
    if (row * column) > 1:
        res_layout = [i for i in range(row * column)]
        # qc_t = transpile(qc, layout_method="dense", coupling_map=cpp)
        # print("done ez layout")
        return res_layout, None
    else:
        qc_t = transpile(
            qc,
            layout_method="sabre",
            coupling_map=cpp,
        )

    print(f"done initial")
    # qc_t = transpile(qc,layout_method="dense",coupling_map=cpp)
    # print(f"by using the sabre layout, we only need {qc_t.count_ops()}")
    for i in qc_t._layout.get_physical_bits():
        temp.append(i)
    return temp, qc_t._layout


def cx_cycle_line_syca(
    linear_list, offset, num_qubit, qaoa_circuit, p2l, G, all_cx_locations, layer
):
    """_summary_

    Args:
        num_qubit (_type_): qubit number
        qc (_type_): quantum circuit
        p2l (_type_): current mapping

    Returns:
        _type_: updated quantum circuit
    """
    for i in range(num_qubit):
        if i % 2 == 0 and i + 1 < num_qubit:
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]

            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_c, Logical_t)] = -layer
                print(
                    f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({linear_list[i]},{linear_list[i+1]})"
                )
                G.remove_edge(Logical_c, Logical_t)
            if G.has_edge(Logical_t, Logical_c):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_t, Logical_c)] = -layer
                print(
                    f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are  ({linear_list[i]},{linear_list[i+1]})"
                )
                G.remove_edge(Logical_t, Logical_c)
    for i in range(num_qubit):
        if i % 2 != 0 and i + 1 < num_qubit:
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]

            if G.has_edge(Logical_c, Logical_t):
                all_cx_locations[(Logical_c, Logical_t)] = -layer
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                print(
                    f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({i},{i+1})"
                )

                G.remove_edge(Logical_c, Logical_t)
            if G.has_edge(Logical_t, Logical_c):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_t, Logical_c)] = -layer
                print(
                    f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({i},{i+1})"
                )

                G.remove_edge(Logical_t, Logical_c)
    return qaoa_circuit


def swap_cycle_line_syca(linear_list, offset, num_qubit, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    for i in range(num_qubit):
        if i % 2 != 0 and i + 1 < num_qubit:
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            print(f"Add swap between {linear_list[i]} and {linear_list[i+1]}")
            update_map(p2l, linear_list[i], linear_list[i + 1])
            # print(f"new l2p {l2p}")

    for i in range(num_qubit):
        if i % 2 == 0 and i + 1 < num_qubit:
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]

            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            print(f"Add swap between {linear_list[i]} and {linear_list[i+1]}")
            update_map(p2l, linear_list[i], linear_list[i + 1])

    return qaoa_circuit


def general_graph_pattern_sycamore(row, column, qc_res, dG, mapping):
    G = nx.Graph(dG)

    ## we are missing initial mapping here,here is the current inital mapping
    physcial_logical = mapping
    # print(physcial_logical)
    ## we later will modified this into a initial mapping

    all_cx_locations = {}
    unit_inter_gate_locations = []
    # qc_res = QuantumCircuit(node_number,node_number)
    # print(type(qc_res))
    # print(f"the graph has {len(G.edges())} edges")
    layer = 0
    ### IE
    unit = 0
    #! row number n*m probelm 15 qubit, 4*4, input should change to n*m, add one more layer,need to fix
    #! add a counter here instead of using the depth of the circuit.
    #### finish  inter A,B
    # print(f"remain {len(G.edges())} edges")
    ##? we here unpdate the row based on the partial
    ##?
    totaln = len(G.edges())
    line_layer = 0
    while unit < row / 2 and len(G.edges()) > 0:
        print(unit)
        save_p2l = physcial_logical.copy()
        save_l2p = {v: k for k, v in save_p2l.items()}
        temp_dic = []
        for edge in G.edges():
            c = save_l2p[edge[0]]
            t = save_l2p[edge[1]]
            c_r = np.floor(c / column)
            t_r = np.floor(t / column)
            print(
                f"physical qubits are{c},{t} in layer {c_r} and {t_r} logical qubits are {edge[0]},{edge[1]}"
            )
            if abs(c_r - t_r) == 1:
                temp_dic.append([c, t])
        unit_inter_gate_locations.append(temp_dic)
        layer = 0

        ####! here we do the even AB CD EF

        while layer < column:
            # print(f'even work {layer}')
            cx_cycle_syca(
                row, column, qc_res, physcial_logical, G, all_cx_locations, layer, unit
            )
            # if len(G.edges()) == totaln - column * column:
            #     break
            swap_cycle_syca(row, column, qc_res, physcial_logical, G, layer)
            layer += 1
        while layer < column * 2:
            # print(f'odd work {layer}')
            ####! here we do the odd  BC DE
            cx_cycle_syca(
                row, column, qc_res, physcial_logical, G, all_cx_locations, layer, unit
            )
            # if len(G.edges()) == totaln - column * column:
            #     break
            swap_cycle_syca(row, column, qc_res, physcial_logical, G, layer)
            layer += 1
        # if unit == row / 2 - 1:
        #     break
        # swap_between_line_syca(row, column, qc_res, physcial_logical, G)
        unit += 1

    offset_line = 0
    #!  we need to deal with intra cx and swap using the linear_pattern
    # line_index = 0
    # linear_list = []
    # for i in range(column):
    #     linear_list.append(i)
    #     linear_list.append(i + column)
    # print(linear_list)
    # if len(G.edges()) > 0:
    #     for i in range(0, row, 2):
    #         # print(f"the offset is {offset_line}")
    #         line_index = 0
    #         while line_index < column:
    #             # print("intra work")
    #             cx_cycle_line_syca(
    #                 linear_list,
    #                 0,
    #                 column * 2,
    #                 qc_res,
    #                 physcial_logical,
    #                 G,
    #                 all_cx_locations,
    #                 line_layer,
    #             )
    #             swap_cycle_line_syca(
    #                 linear_list, 0, column * 2, qc_res, physcial_logical, G
    #             )
    #             line_index += 1
    #             line_layer += 1
    # offset_line += column * 2
    # print(f"remain {len(G.edges())} edges")

    # #print(type(qc_res))
    # print(f"remain {len(G.edges())} edges")

    return qc_res, physcial_logical, all_cx_locations, unit_inter_gate_locations
    # for qaoa_circuit in res_qc_list:
    #     #print(
    #         f"the circuit has depths of {qaoa_circuit.depth()} with gates counts of {qaoa_circuit.count_ops()}"
    #     )
    # for qaoa_circuit in res_qc_list:
    #     #print(
    #         f"the circuit has depths of {qaoa_circuit.depth()} with gates counts of {qaoa_circuit.count_ops()}"
    #     )


def general_graph_pattern_sycamore_partial(
    offsetlist, row, column, qc_res, dG, mapping
):
    G = nx.Graph(dG)
    ##? start row and colum is 0,1

    ##? end row and colum is 2,3

    ## we are missing initial mapping here,here is the current inital mapping
    physcial_logical = mapping
    # print(physcial_logical)
    ## we later will modified this into a initial mapping

    all_cx_locations = {}
    unit_inter_gate_locations = []
    # qc_res = QuantumCircuit(node_number,node_number)
    # print(type(qc_res))
    # print(f"the graph has {len(G.edges())} edges")
    layer = 0
    ### IE
    unit = 0
    #! row number n*m probelm 15 qubit, 4*4, input should change to n*m, add one more layer,need to fix
    #! add a counter here instead of using the depth of the circuit.
    #### finish  inter A,B

    # print(f"remain {len(G.edges())} edges")
    ##? we here unpdate the row based on the partial
    ##?

    line_layer = 0

    while unit < row / 2 and len(G.edges()) > 0:
        # ? our list information is here [[0, 0, 1, 1], [1, 1, 2, 1], None]
        if unit == len(offsetlist):
            break
        cur_info = offsetlist[unit]
        # ? current list information is here [[0, 0, 1, 1], [1, 1, 2, 1], None]
        # save_p2l = physcial_logical.copy()
        # save_l2p = {v: k for k , v in save_p2l.items()}
        # temp_dic = []
        # for edge in G.edges():
        #     c = save_l2p[edge[0]]
        #     t = save_l2p[edge[1]]
        #     c_r = np.floor(c/column)
        #     t_r = np.floor(t/column)
        #     print(f"physical qubits are{c},{t} in layer {c_r} and {t_r} logical qubits are {edge[0]},{edge[1]}")
        #     if  abs(c_r - t_r) == 1:
        #         temp_dic.append([c,t])
        # unit_inter_gate_locations.append(temp_dic)
        layer = 0

        ####! here we do the even AB CD EF

        while layer < column:
            # print(f'even work {layer}')
            ##? for cx gates we do not need to do anything since rest cx gates are not used.
            cx_cycle_syca(
                row, column, qc_res, physcial_logical, G, all_cx_locations, layer, unit
            )
            ##? for swaps we need to modified it using current layer infomations
            # swap_cycle_syca(row,column, qc_res, physcial_logical, G, layer,cur_info)
            swap_cycle_syca_offset(
                row, column, qc_res, physcial_logical, G, layer, cur_info
            )
            layer += 1
        while layer < column * 2:
            # print(f"odd work {layer}")
            ####! here we do the odd  BC DE
            cx_cycle_syca(
                row, column, qc_res, physcial_logical, G, all_cx_locations, layer, unit
            )
            # swap_cycle_syca(row,column, qc_res, physcial_logical, G, layer,cur_info)
            swap_cycle_syca_offset(
                row, column, qc_res, physcial_logical, G, layer, cur_info
            )
            layer += 1

        swap_between_line_syca(row, column, qc_res, physcial_logical, G)
        unit += 1

    offset_line = 0
    #!  we need to deal with intra cx and swap using the linear_pattern
    line_index = 0
    if len(G.edges()) > 0:
        for i in range(0, row, 2):
            # print(f"the offset is {offset_line}")
            line_index = 0
            while line_index < column:
                # print("intra work")
                cx_cycle_line_syca(
                    offset_line,
                    column * 2,
                    qc_res,
                    physcial_logical,
                    G,
                    all_cx_locations,
                    line_layer,
                )
                swap_cycle_line_syca(
                    offset_line, column * 2, qc_res, physcial_logical, G
                )
                line_index += 1
                line_layer += 1
            offset_line += column * 2
        # print(f"remain {len(G.edges())} edges")

        # #print(type(qc_res))
    # print(f"remain {len(G.edges())} edges")

    return qc_res, physcial_logical, all_cx_locations, unit_inter_gate_locations


def cx_cycle_syca(
    row, column, qaoa_circuit, p2l, G, all_cx_locations, layer, unit
) -> QuantumCircuit:
    """_summary_

    Args:
        num_qubit (_type_): qubit number
        qc (_type_): quantum circuit
        p2l (_type_): current mapping

    Returns:
        _type_: updated quantum circuit
    """
    ##### graph informations are logical not physical

    ##### add inter cx gate A,B... C,D
    if layer < column:
        for cur_row in range(0, row - 1, 2):
            for cur_col in range(column):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col

                logical_c = p2l[control_index]
                logical_t = p2l[target_index]
                if G.has_edge(logical_c, logical_t):
                    qaoa_circuit.rzz(0, control_index, target_index)
                    all_cx_locations[(logical_c, logical_t)] = layer + 2 * column * unit
                    # active_gates[(control_index,target_index)] = layer
                    G.remove_edge(logical_c, logical_t)
    ##### add inter cx gate A,B... C,D
    else:
        for cur_row in range(1, row - 1, 2):
            for cur_col in range(column):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col

                logical_c = p2l[control_index]
                logical_t = p2l[target_index]
                if G.has_edge(logical_c, logical_t):
                    qaoa_circuit.rzz(0, control_index, target_index)
                    all_cx_locations[(logical_c, logical_t)] = layer + 2 * column * unit
                    # active_gates[(control_index,target_index)] = layer
                    G.remove_edge(logical_c, logical_t)

    return qaoa_circuit


def swap_cycle_syca(row, column, qaoa_circuit, p2l, G, layer):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    if layer < column:

        ####swap 1 AB,CD we need to deal with partial AB CD
        # ? [[2, 6], [4, 10], None] for AB
        start_index = 0
        for cur_row in range(0, row - 1, 2):

            # ? here is the A,B which is [2,6]
            # column_min = min(cur_list[0]%column,cur_list[1]%column)
            # column_max = max(cur_list[0]%column,cur_list[1]%column)
            # ?here min is 2 max is also 2
            for cur_col in range(layer % 2, column, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap2 should be diagonal
        for cur_row in range(0, row - 1, 2):
            for cur_col in range(1, column):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col - 1
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap3 should be same as swap 1
        for cur_row in range(0, row - 1, 2):
            for cur_col in range(layer % 2, column, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
    else:
        for cur_row in range(1, row - 1, 2):
            for cur_col in range((layer + 1) % 2, column, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap2 should be diagonal
        for cur_row in range(1, row - 1, 2):
            for cur_col in range(0, column - 1):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col + 1
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)

        ####swap3 should be same as swap 1
        for cur_row in range(1, row - 1, 2):
            for cur_col in range((layer + 1) % 2, column, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
    return qaoa_circuit


def swap_cycle_syca_offset(row, column, qaoa_circuit, p2l, G, layer, offset_list):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    if layer < column:

        ####swap 1 AB,CD we need to deal with partial AB CD
        # ? [[0, 0, 1, 1], [1, 1, 2, 1], None] for AB  CD
        start_index = 0
        for cur_row in range(0, row - 1, 2):
            # ? [0, 0, 1, 1]
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the A,B which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            # print(f'the partial colum from{column_min} toi {column_max}')
            # ?here min is 2 max is also 2
            if column_min == column_max:
                # print(f"skip")
                continue
            for cur_col in range(column_min + layer % 2, column_max + 1, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap2 should be diagonal
        for cur_row in range(0, row - 1, 2):
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the A,B which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            # print(f'the partial colum from{column_min} toi {column_max}')
            if column_min == column_max:
                continue
            for cur_col in range(1 + column_min, column_max + 1):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col - 1
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap3 should be same as swap 1
        for cur_row in range(0, row - 1, 2):
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the A,B which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            if column_min == column_max:
                continue
            for cur_col in range(column_min + layer % 2, column_max + 1, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
    else:
        for cur_row in range(1, row - 1, 2):
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the B,C which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            if column_min == column_max:
                continue
            for cur_col in range((layer + 1) % 2 + column_min, column_max + 1, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
        ####swap2 should be diagonal
        for cur_row in range(1, row - 1, 2):
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the A,B which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            if column_min == column_max:
                continue
            for cur_col in range(column_min, column_max):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col + 1
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)

        ####swap3 should be same as swap 1
        for cur_row in range(1, row - 1, 2):
            cur_list = offset_list[cur_row]
            if cur_list == None:
                continue
            # ? here is the A,B which is [2,6]
            column_min = cur_list[1]
            column_max = cur_list[3]
            if column_min == column_max:
                continue
            for cur_col in range((layer + 1) % 2 + column_min, column_max + 1, 2):
                control_index = cur_row * column + cur_col
                target_index = (cur_row + 1) * column + cur_col
                qaoa_circuit.swap(control_index, target_index)
                update_map(p2l, control_index, target_index)
    return qaoa_circuit


def generate_sycamore(n, m):
    res_list = []
    for row in range(n - 1):
        for col in range(m):

            c_index = m * row + col
            t_index = (m) * (row + 1) + col
            temp = [c_index, t_index]
            res_list.append(temp)
    for row in range(1, n):
        if row % 2 == 1:
            for col in range(m - 1):
                c_index = m * row + col
                t_index = (m) * (row - 1) + col + 1
                temp = [c_index, t_index]
                res_list.append(temp)
        else:
            for col in range(1, m):
                c_index = m * row + col
                t_index = (m) * (row - 1) + col - 1
                temp = [c_index, t_index]
                res_list.append(temp)

    return res_list


def read_graph(path):
    res_int = []
    with open(f"{path}") as f:
        lines = f.readlines()
    for l in lines:
        res_int.append(re.findall(r"\b\d+\b", l))

    # print(res_int[0][0])
    GG = nx.Graph()
    num_node = int(res_int[0][0])
    for i in range(num_node):
        GG.add_node(i)
    for i in range(1, len(res_int)):
        l = int(res_int[i][0])
        r = int(res_int[i][1])
        GG.add_edge(l, r)
    return GG, num_node


# res_list = generate_sychmore(4,4)
# dis_g = nx.Graph()
# for i in range(k*k):
#     dis_g.add_node(i)
# for p in res_list:
#     dis_g.add_edge(p[0],p[1])
# print(res_list)
# nx.draw(dis_g,with_labels = True)
# sp = dict(nx.all_pairs_shortest_path_length(dis_g))
# print(sp[0][2])


def get_layout(GG, num):
    """_summary_

    Args:
        GG (_type_): Graph
        num (_type_): Number of nodes

    Returns:
        _type_: layout information using sabre
    """
    qc = QuantumCircuit(num)
    for e in GG.edges():
        qc.cx(e[0], e[1])
    cpp = CouplingMap.from_line(num)
    temp = []
    qc_t = transpile(qc, layout_method="sabre", coupling_map=cpp)

    # qc_t = transpile(qc,layout_method="dense",coupling_map=cpp)
    for i in qc_t._layout.get_physical_bits():
        temp.append(i)
    return temp, qc_t._layout


def analysis_unitary(circuit):
    unitary_circuit = QuantumCircuit(len(circuit.qubits))
    qc_dag = circuit_to_dag(circuit)
    layer = qc_dag.layers()
    cur_index = 0
    total_unitary = 0
    dependency = {}
    pre_layer_cx = []
    cur_layer_cx = []
    unitary_index = {}
    qc = QuantumCircuit(2, name="Cphase_unitary")
    qc.cx(0, 1)
    qc.rz(np.pi / 2, 1)
    qc.cx(1, 0)
    qc.cx(0, 1)
    unitary_gate = qc.to_instruction()

    for i in layer:
        list_node = i["graph"].op_nodes()
        cur_layer_cx = []
        # print(f"start new layer {cur_index}:")
        # print("******")
        # print(list_node)
        for cur_node in list_node:
            if len(cur_node.qargs) > 1:
                control_q = cur_node.qargs[0].index
                target_q = cur_node.qargs[1].index
                # print(f"the {cur_node.name} gate with control qubit {cur_node.qargs[0].index} target qubit {cur_node.qargs[1].index} ")
                # cur_layer_cx.append([control_q,target_q])
                if cur_node.name == "swap":
                    if [control_q, target_q] in pre_layer_cx:
                        total_unitary += 1
                        unitary_index[cur_index - 1] = [control_q, target_q]
                if cur_node.name == "cx":
                    cur_layer_cx.append([control_q, target_q])

        pre_layer_cx = []
        for ccx in cur_layer_cx:

            pre_layer_cx.append(ccx)
        cur_index += 1
    cur_index = 0
    print(unitary_index)
    qc_dag = circuit_to_dag(circuit)
    layer = qc_dag.layers()
    for i in layer:
        list_node = i["graph"].op_nodes()

        # print(f"start new layer {cur_index}:")
        # print("******")
        # print(list_node)
        for cur_node in list_node:
            if len(cur_node.qargs) > 1:
                control_q = cur_node.qargs[0].index
                target_q = cur_node.qargs[1].index
                # print(
                #     f"the {cur_node.name} gate with control qubit {cur_node.qargs[0].index} target qubit {cur_node.qargs[1].index} "
                # )
                # cur_layer_cx.append([control_q,target_q])
                if cur_node.name == "swap":
                    if (
                        cur_index in unitary_index.keys()
                        and [control_q, target_q] == unitary_index[cur_index]
                    ) == False:
                        unitary_circuit.swap(control_q, target_q)
                if cur_node.name == "cx":

                    if (
                        cur_index in unitary_index.keys()
                        and [control_q, target_q] == unitary_index[cur_index]
                    ):
                        print("find!!!")
                        unitary_circuit.append(unitary_gate, [control_q, target_q])
                    else:
                        unitary_circuit.cx(control_q, target_q)

    return total_unitary, unitary_index, unitary_circuit


def run_program_2D(num, GG):
    """_summary_

    Args:
        num (_type_): row/colum num
        GG (_type_): Graph

    Returns:
        _type_: circuit and layout
    """
    # here we use the linear program and the inital mapping is index, we need to update it later, we can update to 2D coupling map
    # cp = get_linear_coup(num)

    cm = CouplingMap.from_grid(num, num)
    print(f"the size is {num}")
    # cm = CouplingMap(cp)
    current_initial_layout, relayout = get_2D_layout(GG, num, cm)
    # print(f"intial logical mapping is {current_initial_layout}")
    # print(f"intial graph is {len(GG.edges())} gates {GG.edges()} ")

    Logical_physical = {}
    for i in range(num * num):
        Logical_physical[i] = current_initial_layout[i]
    physcial_logical = {v: k for k, v in Logical_physical.items()}
    # print(f"correspond phyiscal mapping is {physcial_logical}")

    qaoa_circuit = QuantumCircuit(num * num)
    # print(f"we remain {len(GG.edges())}")
    dG = nx.Graph(GG)
    ## here we need to update our linear pattern to 2D patterns
    (
        qaoa_circuit,
        physcial_logical,
        all_cx_locations,
        res_counts,
    ) = general_graph_pattern_2D(num, qaoa_circuit, dG, physcial_logical)
    print(f"the counts are {res_counts}")
    # print(f"cx locations are{all_cx_locations}")
    # print(qaoa_circuit.draw())
    physcial_logical = {v: k for k, v in Logical_physical.items()}
    # print(f"now mapping is {physcial_logical}")
    # print(f"****start heuristic findings*********")
    # print(f"we remain {len(GG.edges())}")
    res_qc_list = iteration_find(num, GG, physcial_logical, all_cx_locations, cm)
    # for qaoa_circuit in res_qc_list:
    #     #print(
    #         f"the circuit has depths of {qaoa_circuit.depth()} with gates counts of {qaoa_circuit.count_ops()}"
    #     )
    # for qaoa_circuit in res_qc_list:
    #     #print(
    #         f"the circuit has depths of {qaoa_circuit.depth()} with gates counts of {qaoa_circuit.count_ops()}"
    #     )
    return res_qc_list, relayout


def cx_cycle_hex(
    linear_list, outline_list, qaoa_circuit, p2l, G, all_cx_locations, layer, dis_g
) -> QuantumCircuit:
    """_summary_
    s
        Args:
            num_qubit (_type_): qubi number
            qc (_type_): quantum circuit
            p2l (_type_): current mapping

        Returns:
            _type_: updated quantum circuit
    """
    ##### graph informations are logical not physical
    ###! outline_list information need to update
    for k in outline_list:
        temp = [n for n in dis_g.neighbors(k)]
    print(temp)
    for i in range(len(linear_list)):
        if i % 2 != 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({linear_list[i]},{linear_list[i+1]})"
                # )
                G.remove_edge(Logical_c, Logical_t)
    for i in range(len(linear_list)):
        if i % 2 == 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({linear_list[i]},{linear_list[i+1]})"
                # )
                G.remove_edge(Logical_c, Logical_t)
    for k in outline_list:
        #! these are physical qubits
        Logical_c = p2l[k]
        temp = [n for n in dis_g.neighbors(k)]
        for t in temp:
            Logical_t = p2l[t]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, k, t)
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"Outline: we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({k},{t})"
                # )
                G.remove_edge(Logical_c, Logical_t)

    return qaoa_circuit


def cx_cycle_hex_even(
    linear_list, outline_list, qaoa_circuit, p2l, G, all_cx_locations, layer, dis_g
) -> QuantumCircuit:
    """_summary_
    s
        Args:
            num_qubit (_type_): qubi number
            qc (_type_): quantum circuit
            p2l (_type_): current mapping

        Returns:
            _type_: updated quantum circuit
    """
    ##### graph informations are logical not physical
    ###! outline_list information need to update

    for i in range(len(linear_list)):
        if i % 2 == 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({linear_list[i]},{linear_list[i+1]})"
                # )
                G.remove_edge(Logical_c, Logical_t)
    for k in outline_list:
        #! these are physical qubits
        Logical_c = p2l[k]
        temp = [n for n in dis_g.neighbors(k)]
        for t in temp:
            Logical_t = p2l[t]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, k, t)
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"Outline: we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({k},{t})"
                # )
                G.remove_edge(Logical_c, Logical_t)

    return qaoa_circuit


def cx_cycle_hex_out(
    linear_list, outline_list, qaoa_circuit, p2l, G, all_cx_locations, layer, dis_g
) -> QuantumCircuit:
    """_summary_
    s
        Args:
            num_qubit (_type_): qubi number
            qc (_type_): quantum circuit
            p2l (_type_): current mapping

        Returns:
            _type_: updated quantum circuit
    """
    ##### graph informations are logical not physical
    ###! outline_list information need to update
    for k in outline_list:
        temp = [n for n in dis_g.neighbors(k)]
    for k in outline_list:
        #! these are physical qubits
        Logical_c = p2l[k]
        temp = [n for n in dis_g.neighbors(k)]
        for t in temp:
            Logical_t = p2l[t]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, k, t)
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"Outline: we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({k},{t})"
                # )
                G.remove_edge(Logical_c, Logical_t)

    return qaoa_circuit


def cx_cycle_hex_odd(
    linear_list, outline_list, qaoa_circuit, p2l, G, all_cx_locations, layer, dis_g
) -> QuantumCircuit:
    """_summary_
    s
        Args:
            num_qubit (_type_): qubi number
            qc (_type_): quantum circuit
            p2l (_type_): current mapping

        Returns:
            _type_: updated quantum circuit
    """
    ##### graph informations are logical not physical
    ###! outline_list information need to update
    for i in range(len(linear_list)):
        if i % 2 != 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            if G.has_edge(Logical_c, Logical_t):
                qaoa_circuit.rzz(0, linear_list[i], linear_list[i + 1])
                all_cx_locations[(Logical_c, Logical_t)] = layer
                # print(
                #     f"we finsh logical cx({Logical_c},{Logical_t}) which in physical are ({linear_list[i]},{linear_list[i+1]})"
                # )
                G.remove_edge(Logical_c, Logical_t)
    return qaoa_circuit


def swap_cycle_hex_odd(linear_list, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    for i in range(len(linear_list)):
        if i % 2 != 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            # print(f"original l2p {l2p}")
            update_map(p2l, linear_list[i], linear_list[i + 1])

    return qaoa_circuit


def swap_cycle_hex_even(linear_list, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return

    for i in range(len(linear_list)):
        if i % 2 == 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            # print(f"original l2p {l2p}")
            update_map(p2l, linear_list[i], linear_list[i + 1])
    return qaoa_circuit


def heuristic(qc_res, physcial_logical, G):
    for e in G.edges():
        qc_res.cx(e[0], e[1])
        G.remove_edge(e[0], e[1])
    return qc_res


def swap_cycle_hex_1(linear_list, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    for i in range(len(linear_list)):
        if i % 2 == 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            # print(f"original l2p {l2p}")
            update_map(p2l, linear_list[i], linear_list[i + 1])
    return qaoa_circuit


def general_graph_pattern_hex(
    row, column, line_list, outline_list, qc_res, dG, mapping, dis_g
):
    G = nx.Graph(dG)
    ## we are missing initial mapping here,here is the current inital mapping
    physcial_logical = mapping
    # print(physcial_logical)
    ## we later will modified this into a initial mapping
    all_cx_locations = {}
    # qc_res = QuantumCircuit(node_number,node_number)
    # print(type(qc_res))
    # print(f"the graph has {len(G.edges())} edges")
    layer = 0
    #! we first finish the cx gates for line
    # print(f"remain {len(G.edges())} edges")
    ##? we here unpdate the row based on the partial
    ##?

    for i in range(int(np.ceil(len(line_list) / 2 + 1))):
        # print(f"*******layer{layer}")
        cx_cycle_hex_out(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        cx_cycle_hex_odd(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )

        swap_cycle_hex_odd(line_list, qc_res, physcial_logical, G)
        cx_cycle_hex_out(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        cx_cycle_hex_even(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        swap_cycle_hex_even(line_list, qc_res, physcial_logical, G)
        layer += 1
        # if layer > 0:
        #     return qc_res, physcial_logical, all_cx_locations
    ##? here we do the swap for outline
    for k in outline_list:

        temp = [n for n in dis_g.neighbors(k)]
        temp.sort()
        kk = temp[0]
        qc_res.swap(k, kk)
        # print(f"original l2p {l2p}")
        update_map(physcial_logical, k, kk)
    ##? then we do the linear again
    for i in range(int(np.ceil(len(line_list) / 2) + 2)):
        # print(f"*******layer{layer}")
        if len(G.edges()) == 0:
            break
        cx_cycle_hex_out(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        cx_cycle_hex_odd(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )

        swap_cycle_hex_odd(line_list, qc_res, physcial_logical, G)
        cx_cycle_hex_out(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        cx_cycle_hex_even(
            line_list,
            outline_list,
            qc_res,
            physcial_logical,
            G,
            all_cx_locations,
            layer,
            dis_g,
        )
        swap_cycle_hex_even(line_list, qc_res, physcial_logical, G)
        # cx_cycle_hex(line_list,outline_list,qc_res,physcial_logical,G,all_cx_locations,layer)
        # swap_cycle_hex_1(line_list,qc_res,physcial_logical,G)
        layer += 1

    heuristic(qc_res, physcial_logical, G)
    # print(f"remain {len(G.edges())} edges")
    # for i in range(0,row,2):
    #     print(f"the offset is {offset_line}")

    #     while line_index < column:
    #         print("intra work")
    #         cx_cycle_line_syca(offset_line,column*2,qc_res,physcial_logical,G,all_cx_locations,line_index)
    #         swap_cycle_line_syca(offset_line,column*2,qc_res,physcial_logical,G)
    #         line_index += 1

    #     offset_line += column*2

    # #print(type(qc_res))

    return qc_res, physcial_logical, all_cx_locations


def swap_cycle_hex_0(linear_list, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    for i in range(len(linear_list)):
        if i % 2 != 0 and i + 1 < len(linear_list):
            Logical_c = p2l[linear_list[i]]
            Logical_t = p2l[linear_list[i + 1]]
            qaoa_circuit.swap(linear_list[i], linear_list[i + 1])
            # print(f"original l2p {l2p}")
            update_map(p2l, linear_list[i], linear_list[i + 1])
    return qaoa_circuit


def generate_linear_list(row, column, ro, co):
    linear_list = []
    x = co
    reverse = True
    for r in range(row):
        if r % 2 == 0:
            left = int(r / 2) * (column + x + 1)
            right = int(r / 2) * (column + x + 1) + column

            temp_list = [i for i in range(left, right)]
            # print(temp_list)
            if reverse:
                temp_list.reverse()
                for i in temp_list:
                    linear_list.append(i)
                reverse = False
            else:
                for i in temp_list:
                    linear_list.append(i)
                reverse = True
        else:
            if int((r - 1) / 2) == 0:

                linear_list.append(int((r + 1) / 2) * (column + x + 1) - x - 1)
            else:
                linear_list.append(int((r + 1) / 2) * (column + x + 1) - 1)
    return linear_list


def create_qaoa_circ(G, theta):

    """
    Creates a parametrized qaoa circuit

    Args:
        G: networkx graph
        theta: list
               unitary parameters

    Returns:
        qc: qiskit circuit
    """

    nqubits = len(G.nodes())
    p = len(theta) // 2  # number of alternating unitaries
    qc = QuantumCircuit(nqubits)

    beta = theta[:p]
    gamma = theta[p:]

    # initial_state
    for i in range(0, nqubits):
        qc.h(i)

    for irep in range(0, p):

        # problem unitary
        for pair in list(G.edges()):
            qc.rzz(2 * gamma[irep], pair[0], pair[1])

        # mixer unitary
        for i in range(0, nqubits):
            qc.rx(2 * beta[irep], i)

    qc.measure_all()

    return qc


def generate_hh(row, column, ro, co):
    res_list = []
    k = co
    res_list = []
    for r in range(row):

        if r % 2 == 0:
            for c in range(column - 1):

                c_i = int((r / 2) * (column + k + 1) + c)
                t_i = c_i + 1
                temp = [c_i, t_i]

                res_list.append(temp)
        else:
            for c in range(k + 1):
                c_i = int((r - 1) / 2 * (column + k + 1) + column + c)
                t_i = int((r - 1) / 2 * (column + k + 1) + 4 * (c))
                t_ii = int((r + 1) / 2 * (column + k + 1) + 4 * (c))
                temp = [c_i, t_i]
                res_list.append(temp)
                temp = [c_i, t_ii]
                res_list.append(temp)
    return res_list


def get_size(row, column, ro, co):
    if row % 2 == 0:
        return int((row + 1) / 2) * (column + co + 1)
    else:
        return int((row + 1) / 2) * (column + co + 1) - co - 1


def swap_between_line_syca(row, column, qaoa_circuit, p2l, G):
    """_summary_

    Args:
        num_qubit (_type_): _description_
        qc (_type_): _description_
        p2l (_type_): _description_
        l2p (_type_): _description_
        G (_type_): _description_

    Returns:
        _type_: _description_
    """
    if len(G.edges) == 0:
        return
    ## layer 1
    for cur_row in range(0, row - 1, 2):
        for cur_col in range(column):
            control_index = cur_row * column + cur_col
            target_index = (cur_row + 1) * column + cur_col
            qaoa_circuit.swap(control_index, target_index)
            update_map(p2l, control_index, target_index)
    ## layer 2
    for cur_row in range(1, row - 1, 2):
        for cur_col in range(column):
            control_index = cur_row * column + cur_col
            target_index = (cur_row + 1) * column + cur_col
            qaoa_circuit.swap(control_index, target_index)
            update_map(p2l, control_index, target_index)

    return qaoa_circuit


def run_program_hh(nrow):
    GG = nx.complete_graph(nrow * 2)
    num = 2 * nrow
    # row = 2
    # column = nrow
    # GG, num = read_graph(f"{cur_path}")
    # ro = 1
    # co = 1
    # column = 4 * co + 3
    # row = 2 * ro + 1
    # num = get_size(row, column, ro, co)
    # print(num)
    flag = True
    ro = 1
    co = 1
    column = 4 * co + 3
    row = 2 * ro + 1
    num = get_size(row, column, ro, co)
    while num < len(GG.nodes):
        if flag:
            ro += 1
            row = 2 * ro + 1
            column = 4 * co + 3
            flag = False
        else:
            co += 1
            row = 2 * ro + 1
            column = 4 * co + 3
            flag = True
        num = get_size(row, column, ro, co)
    # print(num)
    res_list = generate_hh(row, column, ro, co)
    print(res_list)

    dis_g = nx.Graph()
    for i in range(num):
        dis_g.add_node(i)
    for p in res_list:
        dis_g.add_edge(p[0], p[1])

    nx.draw(dis_g, with_labels=True)
    sp = dict(nx.all_pairs_shortest_path_length(dis_g))

    cm = CouplingMap(res_list)
    # current_initial_layout, relayout = get_2D_layout(GG, row, column, ro, co, cm)
    qcc = create_qaoa_circ(GG, [1.0, 1.0])
    # res_circuit = transpile(res_circuit, backend_qaoa, optimization_level=3)

    # hmapper = HeuristicMapper(
    #     qcc.qasm(), coupling_map=res_list
    # )  # The mapping pass in 2QAN
    # # qcc_t =transpile(qcc,backend,)
    # init_map = hmapper.run_qiskit(max_iterations=5)
    # # current_initial_layout = [i for i in range(num)]
    linear_list = generate_linear_list(row, column, ro, co)
    outline_list = []
    for i in range(num):
        if i not in linear_list:
            outline_list.append(i)
    # rint(f"the outline qubits are {outline_list}")
    # print(linear_list)
    init_map = [i for i in range(num)]
    Logical_physical = {}

    for i in range(num):
        Logical_physical[i] = init_map[i]
    # for i in range(num):
    #      Logical_physical[i] = i
    physcial_logical = {v: k for k, v in Logical_physical.items()}
    # print(f"correspond phyiscal mapping is {physcial_logical}")
    qaoa_circuit = QuantumCircuit(num)

    qaoa_circuit, p2l, all_cx_locations = general_graph_pattern_hex(
        row,
        column,
        linear_list,
        outline_list,
        qaoa_circuit,
        GG,
        physcial_logical,
        dis_g,
    )
    # print(
    #         f"the best circuit has depths of {qaoa_circuit.depth()} with gates counts of {qaoa_circuit.count_ops()} has depth{qaoa_circuit.depth()}"
    #     )
    return qaoa_circuit


def run_program(nrow):
    # GG, num = read_graph(cur_path)
    print("start pattern generate")
    GG = nx.complete_graph(nrow * 2)
    num = 2 * nrow
    row = 2
    column = nrow
    res_list = generate_sycamore(row, column)
    # dis_g = nx.Graph()
    # for i in range(row * column):
    #     dis_g.add_node(i)
    # for p in res_list:
    #     dis_g.add_edge(p[0], p[1])
    # sp = dict(nx.all_pairs_shortest_path_length(dis_g))
    print(res_list)
    cm = CouplingMap(res_list)
    # current_initial_layout, relayout = get_2D_layout(GG, row, column, cm)
    current_initial_layout = [i for i in range(row * column)]
    Logical_physical = {}
    num = row * column
    for i in range(num):
        Logical_physical[i] = current_initial_layout[i]
    physcial_logical = {v: k for k, v in Logical_physical.items()}
    # print(f"correspond phyiscal mapping is {physcial_logical}")
    qaoa_circuit = QuantumCircuit(num)
    dG = nx.complete_graph(num)
    (
        qc,
        ef,
        all_cx_locations,
        unit_inter_gate_locations,
    ) = general_graph_pattern_sycamore(row, column, qaoa_circuit, GG, physcial_logical)

    print("we use this")
    return qc
    # GG, num = read_graph(cur_path)
    # num = row * column

    # ll = 0
    # offset_list = []
    # for i in unit_inter_gate_locations:
    #     temp_list = [None] * (column)
    #     # print(f"in unit {ll}, involved gates are {i}")
    #     for g in i:
    #         c = g[0]
    #         t = g[1]
    #         row_lar = max(int(np.floor(c / column)), int(np.floor(t / column)))
    #         row_smm = min(int(np.floor(c / column)), int(np.floor(t / column)))
    #         col_lar = max(c % column, t % column)
    #         col_smm = min(c % column, t % column)
    #         if temp_list[row_smm] == None:
    #             temp_list[row_smm] = [row_smm, col_smm, row_lar, col_lar]
    #         else:

    #             temp_list[row_smm] = [
    #                 min(row_smm, temp_list[row_smm][0]),
    #                 min(col_smm, temp_list[row_smm][1]),
    #                 max(row_lar, temp_list[row_smm][2]),
    #                 max(col_lar, temp_list[row_smm][3]),
    #             ]
    #         # print(temp_list)
    #         # print(f'the gate ({c},{t}) are in layer {np.floor(c/m)} and layer {np.floor(t/m)}')
    #     ll += 1
    #     offset_list.append(temp_list)
    # # ? so the offset should be a nested list first [unit0,unit1] then it will be different unit layer
    # # print(unit_inter_gate_locations)
    # # print(offset_list)
    # # Logical_physical = {}

    # for i in range(num):
    #     Logical_physical[i] = current_initial_layout[i]
    # physcial_logical = {v: k for k, v in Logical_physical.items()}
    # # print(f"correspond phyiscal mapping is {physcial_logical}")
    # qaoa_circuit = QuantumCircuit(num)
    # (
    #     qc,
    #     ef,
    #     all_cx_locations,
    #     unit_inter_gate_locations,
    # ) = general_graph_pattern_sycamore_partial(
    #     offset_list, row, column, qaoa_circuit, GG, physcial_logical
    # )

    # print(f"we need in total {qaoa_circuit.depth()} layers")


def get_best_circuit(circuit):

    score = circuit[0].depth()
    res_c = circuit[0]
    for c in circuit:
        cur_score = c.depth()
        if cur_score < score:
            score = cur_score
            res_c = c
    return res_c


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path",
        type=str,
        default="benchmarks/twolocale",
        help="path of the circuit",
    )
    parser.add_argument(
        "--check",
        type=bool,
        default=False,
        help="check correctness",
    )
    parser.add_argument(
        "--number",
        type=int,
        default=10,
        help="2 * n pattern",
    )
    parser.add_argument(
        "--output",
        type=bool,
        default=False,
        help="Output Schedule",
    )
    parser.add_argument(
        "--arch",
        type=str,
        default="syca",
        help="architecture design",
    )
    parser.add_argument(
        "--size",
        type=str,
        default="2,5",
        help="size",
    )
    res_sum = {}
    # n=5
    # qc = run_program(n)

    args = parser.parse_args()
    path_str = args.path
    check = args.check
    output = args.output
    arch = args.arch
    row_number = args.number
    edges_n = row_number * row_number
    if arch == "syca":
        qc = run_program(row_number)
    elif arch == "heavyhex":
        edges_n = row_number * (2 * row_number - 1)
        qc = run_program_hh(row_number)
    else:
        print("error, we do not have this machine")
        return
    # qc = run_program(row_number)
    # qc2 = run_program_hh(row_number)
    if check:
        print("start check the correctness")
        print(qc.count_ops())
        print(qc.depth())

        print(edges_n)
        counts = qc.count_ops()["rzz"]
        if counts != edges_n:
            print("We find errors, the rzz gates not match the edges of graph")
        else:
            print(f"we finish generating {arch} architecture")
            print(
                f"There are in total {edges_n} edges of graph, which matches {counts} rzz gates"
            )
    if output:
        print("Output the circuit scheduling")
        qc_qasm = qc.qasm()
        with open("output.qasm", "w+") as f:
            f.writelines(qc_qasm)
    # print(qc.count_ops())

    # print(row_number)


if __name__ == "__main__":
    main()
