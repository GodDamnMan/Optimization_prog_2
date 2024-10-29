import numpy
from numpy.linalg import norm
from task1 import Simplex
from itertools import *
import random
import math

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

    

    def find_initial_solution(self):
        initial_solution = [1 for _ in range(self.n)]

        for i in range(len(self.constraints)):
            initial_solution.append(self.right_hand_side[i] - sum(self.constraints[i][:self.n]))


        dirty_bit = True
        while dirty_bit:
            dirty_bit = False
            for i in range(len(initial_solution)):
                if initial_solution[i] < 0:
                    dirty_bit = True
                    if i < self.n:
                        new_initial_solution = [(initial_solution[j] if j!=i else 0) for j in range(self.n)]
                        
                        for j in range(len(self.constraints)):
                            a = self.constraints[j][:self.n]
                            b = [a[k] * new_initial_solution[k] for k in range(len(a))]
                            new_initial_solution.append(self.right_hand_side[j] - sum(b))
                        initial_solution = new_initial_solution
                        break
                    
                    
                    ind = i - self.n
                    c = -self.constraints[ind][i]
                    grad = [self.constraints[ind][j]/c if j != i else 0  for j in range(len(self.constraints[ind]))] # + [self.right_hand_side[ind]]
                    s = (sum([i**2 for i in grad])) ** 0.5
                    step = 0.1 - initial_solution[i]
                    grad = [j * step / s  for j in grad]
                    new_initial_solution = [initial_solution[j] + grad[j] for j in range(self.n)]

                    for j in range(len(self.constraints)):
                        a = self.constraints[j][:self.n]
                        b = [a[k] * new_initial_solution[k] for k in range(len(a))]
                        new_initial_solution.append(self.right_hand_side[j] - sum(b))
                    initial_solution = new_initial_solution
                    break

        self.solution_l1 = initial_solution.copy()
        self.solution_l2 = initial_solution.copy()

    def iterior_point_method(self):
        self.find_initial_solution()
        try:
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
        except ValueError:
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
            raise ValueError
        
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
    degeneracy_flag = False
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
        degeneracy_flag = True
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
        if not degeneracy_flag: print("ERROR: NOT A NUMBER")
        else: 
            print('DEGENERACY CASE')
            
        return
    

    except numpy.linalg.LinAlgError:
        print("The method is not applicable!")
        return


userInput()


