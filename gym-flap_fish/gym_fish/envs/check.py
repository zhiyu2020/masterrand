# import test_2d_flow_stream_around_fish_test_pybind as test_fish
import test_2d_flow_stream_around_fish_flap_pybind as train
import test_2d_flow_stream_around_fish_flap_test_pybind as test

r = train.from_sph_relaxation(0)
a = train.from_sph_reload_and_train(1)
b = test.from_sph_reload_and_test(101)

a.SetFreq(3)
a.SetLambda(5)
for i in range(32):
	pressure = a.GetFishPressurePoint(i)
	v_x = a.GetFishVelocityX(i)
	v_y = a.GetFishVelocityY(i)
	pos_x = a.GetFishPositionX(i)
	pos_y = a.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
a.RunCase(1, 0.2)
for i in range(32):
	pressure = a.GetFishPressurePoint(i)
	v_x = a.GetFishVelocityX(i)
	v_y = a.GetFishVelocityY(i)
	pos_x = a.GetFishPositionX(i)
	pos_y = a.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
for i in range(32):
	pressure = b.GetFishPressurePoint(i)
	v_x = b.GetFishVelocityX(i)
	v_y = b.GetFishVelocityY(i)
	pos_x = b.GetFishPositionX(i)
	pos_y = b.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
b.RunCase(101, 0.2)
for i in range(32):
	pressure = b.GetFishPressurePoint(i)
	v_x = b.GetFishVelocityX(i)
	v_y = b.GetFishVelocityY(i)
	pos_x = b.GetFishPositionX(i)
	pos_y = b.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)

for i in range(32):
	pressure = a.GetFishPressurePoint(i)
	v_x = a.GetFishVelocityX(i)
	v_y = a.GetFishVelocityY(i)
	pos_x = a.GetFishPositionX(i)
	pos_y = a.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
a.RunCase(1, 0.4)
for i in range(32):
	pressure = a.GetFishPressurePoint(i)
	v_x = a.GetFishVelocityX(i)
	v_y = a.GetFishVelocityY(i)
	pos_x = a.GetFishPositionX(i)
	pos_y = a.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
b.RunCase(101, 0.3)



