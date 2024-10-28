import numpy
from numpy.linalg import norm
from task1 import Simplex
from itertools import*

class IteriorPoint:

    def __init__(self, coefs:list, constraints:list, right_hand_side:list, epsilon:int, n:int, isMax:bool):
        self.solution_l1 = None
        self.solution_l2 = None
        self.epsilon = epsilon
        self.right_hand_side = right_hand_side
        self.constraints = constraints
        self.coefs = coefs
        self.alpha1 = 0.5
        self.alpha2 = 0.9
        self.iter_val = 1
        self.n = n
        self.isMax = isMax

    def find_all_axis(self):
        all_axis = []
        axe = []
        for i in range(self.n-1):
            axe.append(0)
        axe.append(1)
        all_axis.append(axe)
        for i in range(self.n-1):
            tmp_axe = all_axis[i]
            rotated_axe = tmp_axe[1:] + tmp_axe[:1]
            all_axis.append(rotated_axe)
        return all_axis

    # TODO оно не работает для размерности больше чем 2
    def centroid(self, vertices):
        x, y = 0, 0
        n = len(vertices)
        signed_area = 0
        for i in range(len(vertices)):
            x0, y0 = vertices[i]
            x1, y1 = vertices[(i + 1) % n]
            area = (x0 * y1) - (x1 * y0)
            signed_area += area
            x += (x0 + x1) * area
            y += (y0 + y1) * area
        signed_area *= 0.5
        x /= 6 * signed_area
        y /= 6 * signed_area
        return x, y


    def intercection(self, intercection_try):
        tmp_constraints = []
        for constraint_line in intercection_try:
            tmp_list = []
            for constraint_element in range(self.n):
                tmp_list.append(constraint_line[constraint_element])
            tmp_constraints.append(tmp_list)
        tmp_rhs = self.right_hand_side
        tmp_rhs.append(0)

        intercection_points = []
        for i in range(len(tmp_constraints)-1):
            for j in range(i+1, len(tmp_constraints)):
                M = numpy.matrix([tmp_constraints[i], tmp_constraints[j]])
                C = numpy.matrix([tmp_rhs[i], tmp_rhs[j]])
                el = M.I.dot(C.T).tolist()
                dot = el[0] +  el[1]
                intercection_points.append(list(dot))
        return intercection_points

    def find_initial_solution(self):
        initial_solution = []

        flag = False
        for i in self.right_hand_side:
            if(i < 0):
                flag = True
                break

        if(not flag):
            for _ in range(self.n):
                initial_solution.append(0.00001)
        
            for i in range(len(self.constraints)):
                initial_solution.append(self.right_hand_side[i] - 0.00001 * sum(self.constraints[i][:self.n]))
        
        else:
            all_axis = self.find_all_axis()
            for i in all_axis:
                intersection_try = self.constraints + [i]
                intercection_points = self.intercection(intersection_try)
                # TODO оно не работает для размерности больше чем 2 (centroid)
                x1, x2 = self.centroid(intercection_points)
                flag = True
                for i in range(len(self.constraints)):
                    res_to_check = self.constraints[i][0] * x1 + self.constraints[i][1] * x2
                    if(res_to_check > self.right_hand_side[i]):
                        flag = False
                        break
                    
                if(flag):
                    initial_solution.append(x1)
                    initial_solution.append(x2)
        #TODO self.constraints стоит опасно может улететь
                    for i in range(len(self.constraints)):
                        res = self.right_hand_side[i]
                        for j in range(self.n):
                            res -= initial_solution[j] * self.constraints[i][j]
                        initial_solution.append(res)
        
        print(initial_solution)
        self.solution_l1 = initial_solution.copy()
        self.solution_l2 = initial_solution.copy()

    def iterior_point_method(self):
        self.find_initial_solution()
        while True:
            #calculating previous solution and diagonal matrix for lambda1 and lambda2
            prev_sol_l1 = self.solution_l1
            diagonal_l1 = numpy.diag(self.solution_l1)    

            prev_sol_l2 = self.solution_l2
            diagonal_l2 = numpy.diag(self.solution_l2)

            #calculating matrices for lambda1 and lambda2
            AA_l1 = numpy.dot(self.constraints,diagonal_l1)
            cc_l1 = numpy.dot(diagonal_l1, self.coefs)

            AA_l2 = numpy.dot(self.constraints,diagonal_l2)
            cc_l2 = numpy.dot(diagonal_l2, self.coefs)


            I = numpy.eye(len(self.coefs))


            #components to calculate P for lambda1 and lambda2
            F_l1 = numpy.dot(AA_l1, numpy.transpose(AA_l1))
            FI_l1 = numpy.linalg.inv(F_l1)
            H_l1 = numpy.dot(numpy.transpose(AA_l1), FI_l1)
            
            
            F_l2 = numpy.dot(AA_l2, numpy.transpose(AA_l2))
            FI_l2 = numpy.linalg.inv(F_l2)
            H_l2 = numpy.dot(numpy.transpose(AA_l2), FI_l2)

            #calculating P for lambda1 and lambda2
            P_l1 = numpy.subtract(I, numpy.dot(H_l1, AA_l1))
            P_l2 = numpy.subtract(I, numpy.dot(H_l2, AA_l2))

            #calculating cp for lambda1 and lambda2
            cp_l1 = numpy.dot (P_l1, cc_l1)
            cp_l2 = numpy.dot (P_l2, cc_l2)


            
            nu_l1 = numpy.absolute(numpy.min(cp_l1))
            y_l1 = numpy.add(numpy.ones(len(self.coefs), float), (self.alpha1/nu_l1) * cp_l1)

            nu_l2 = numpy.absolute(numpy.min(cp_l2))
            y_l2 = numpy.add(numpy.ones(len(self.coefs), float), (self.alpha2/nu_l2) * cp_l2)
            


            yy_l1 = numpy.dot(diagonal_l1, y_l1)
            yy_l2 = numpy.dot(diagonal_l2, y_l2)
            
            self.solution_l1 = yy_l1
            self.solution_l2 = yy_l2
            
            if self.iter_val==1 or self.iter_val == 2 or self.iter_val == 3 or self.iter_val == 4:
                x_l1 = "[ "
                x_l2 = "[ "
                for item in self.solution_l1:
                    x_l1 += f'{item:.{self.epsilon}f}' + " "
                x_l1 += "]"
                for item in self.solution_l2:
                    x_l2 += f'{item:.{self.epsilon}f}' + " "
                x_l2 += "]"

                print("In iteration ", self.iter_val, "\n\tfor lambda = 0.5 we have x = ", x_l1, "\n\tfor lambda = 0.9 we have x = ", x_l2)
                self.iter_val = self.iter_val + 1
            
            if norm(numpy.subtract(yy_l1, prev_sol_l1), ord = 2)< 0.00001 or norm(numpy.subtract(yy_l2, prev_sol_l2), ord = 2)< 0.00001:
                break
        x_l1 = "[ "
        x_l2 = "[ "
        for item in self.solution_l1:
            x_l1 += f'{item:.{self.epsilon}f}' + " "
        x_l1 += "]"
        for item in self.solution_l2:
            x_l2 += f'{item:.{self.epsilon}f}' + " "
        x_l2 += "]"
        print ("In the last iteration ", self.iter_val, "\n\tfor lambda = 0.5 we have x = " , x_l1, "\n\tfor lambda = 0.9 we have x = ", x_l2)
        res_l1 = 0
        res_l2 = 0
        for i in range(len(self.coefs)):
            res_l1 += self.coefs[i] * self.solution_l1[i]
            res_l2 += self.coefs[i] * self.solution_l2[i]
        print("\nOptimum:\n\tfor lambda = 0.5 we have:", round(res_l1, self.epsilon), "\n\tfor lambda = 0.9 we have", round(res_l2, self.epsilon))

    def solve_max(self):
        self.iterior_point_method()

    def solve_min(self):
        for item in self.coefs:
            item *= -1
        self.iterior_point_method()

##
# function that will get data from user and check it for some typos, etc.
# and will envoke 2 methods (simplex, iterior-point) if exerything is fine =)
##
def userInput():
    type = input("Greetings, this programm will solve your LP problem using interior-point method.\nEnter the type of the problem(Max/Min): ").lower()
    m = type.index("m")
    type = type[m:m+3]
    if (type != "max" and type != "min"):
        print("ERROR: UNKNOWN TYPE")
        return
    print("NOTE: in input there shouldn't be any slack variables")
    try:
        
        objective_function = input("Enter the coefficients of the objective function: ").split(" ")
        for i in range(objective_function.count("")):
            objective_function.remove("")

        objective_function = list(map(float, objective_function))
        if(len(objective_function) == 0):
            print("ERROR: NO COEFFICIENTS")
            return
        n = len(objective_function)
        

        amount = int(input("Enter amount of the constraints(not assuming x>=0): "))
        if(amount < 1):
            print("ERROR: AMOUNT < 1 ?!")
        constraints = []
        for i in range(amount):
            constraint = input(f"Enter the {i+1} constraint function coefficients: ").split(" ")
            for i in range(constraint.count("")):
                constraint.remove("")
            constraint = list(map(float, constraint))
            if(len(constraint) == 0):
                print("ERROR: NOT ENOUGH COEFFICIENTS")
                return
            constraints.append(constraint)
        

        # Unbounded check
        for j in range(len(constraints[0])):
            f = True
            for i in range(len(constraints)):
                if constraints[i][j] > 0:
                    f = False
                    continue
            if f:
                print("The problem does not have solution!")
                return

        right_hand_side = input("Enter the right-hand side numbers: ").split(" ")
        for i in range(right_hand_side.count("")):
            right_hand_side.remove("")
        right_hand_side = list(map(float, right_hand_side))
        flag_rhs = False
        if(len(right_hand_side) != amount):
            print("ERROR: NOT ENOUGH COEFFICIENTS")
            return
        else:
            for i in right_hand_side:
                if(i < 0):
                    print("The  Simplex method is not applicable!(Input contain negative RHS value)")
                    flag_rhs = True

        accuracy = input("Enter the approximation accuracy: ")

        if accuracy.count('.') == 0: 
            accuracy = int(accuracy)
        else:
            k = 0
            accuracy = float(accuracy)
            while True:
                k += 1
                accuracy *= 10
                if accuracy > 0:
                    break
            accuracy = k

        if(not flag_rhs):
            lp_simplex = Simplex(objective_function, constraints, right_hand_side, accuracy)

        #the following block of code make change of input data for different methods
        for i in range(len(constraints)):
            objective_function.append(0)
            constraint_iter = constraints[i]
            for j in range(len(constraints)):
                if i == j:
                    constraint_iter.append(1)
                else:
                    constraint_iter.append(0)

        # lp_iterior_point = IteriorPoint(objective_function, constraints, right_hand_side, accuracy, n)
    
        if(type == "max"):
            lp_iterior_point = IteriorPoint(objective_function, constraints, right_hand_side, accuracy, n, True)
            if(not flag_rhs):
                print("===============================================================")
                print ("SIMPLEX")
                lp_simplex.solve_maximize()
            print("===============================================================")
            print ("INTERIOR POINT")
            lp_iterior_point.solve_max()
            print("===============================================================")


        else:
            lp_iterior_point = IteriorPoint(objective_function, constraints, right_hand_side, accuracy, n, False)
            if(not flag_rhs):
                print("===============================================================")
                print ("SIMPLEX")
                lp_simplex.solve_minimize()
            print("===============================================================")
            print ("INTERIOR POINT")
            lp_iterior_point.solve_min()
            print("===============================================================")
        

    except ValueError:
        print("ERROR: NOT A NUMBER")
        return
#TODO валится на degeneracy alternative optima что то проходит что то валится но вроде так и должно быть
    except numpy.linalg.LinAlgError:
        print("The method is not applicable!")
        return



userInput()
