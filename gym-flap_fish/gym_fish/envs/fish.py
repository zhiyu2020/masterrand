import math
import re
import time
import torch
import io
import numpy as np
import test_2d_flow_stream_around_fish_flap_pybind as fish
import test_2d_flow_stream_around_fish_flap_test_pybind as test
import gymnasium as gym
from gymnasium import spaces

def checknan(reward):
    if math.isnan(reward):
        return True
    return False


class FISHEnv(gym.Env):
    metadata = {"render_modes": ["human", "rgb_array"], "render_fps": 30}

    def __init__(self, render_mode=None):
        self.tank_lenth = 0.64
        self.tank_height = 0.4
        self.alpha = 0.1
        self.episode = 1
        self.output_episode = 100
        self.time_per_action = 0.1
        # self.goal_x = 0.3 # back of the flap
        # self.goal_y = 0.45 # back of the flap
        self.goal_x = 0.1 * 0.5# back of the flap
        self.goal_y = 0.3 # back of the flap
        low_action = np.array([-1, -1]).astype(np.float32)
        high_action = np.array([1, 1]).astype(np.float32)
        # 32 * 5 pressurepoint, V_x, V_y, pos_x, pos_y
        low_obs = np.full(32, -10).astype(np.float32)
        high_obs = np.full(32, 10).astype(np.float32)
        self.obs = np.zeros(32)
        # arr ????
        self.arr = np.full(32, 100)

        self.action_space = spaces.Box(low_action, high_action)
        self.observation_space = spaces.Box(low_obs, high_obs)

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)

        self.total_reward = 0.0
        self.action_time = 0.0
        self.action_time_steps = 0
        self.frequency = 3
        self.lamb = 2

        self.train_fish = fish.from_sph_relaxation(self.episode) # generate configuration xml fish file

        self.train_fish = fish.from_sph_reload_and_train(self.episode)
        self.train_fish.SetFreq(self.frequency)
        self.train_fish.SetLambda(self.lamb)
        self.train_fish.RunCase(self.episode, self.action_time)
        # total 32 points for corresponding point
        self.fish_start_x = self.train_fish.GetFishPositionX(0)
        self.fish_start_y = self.train_fish.GetFishPositionY(0)

        for i in range(32):
            self.obs[i] = self.train_fish.GetFishPressurePoint(i) / self.arr[i]
        # self.obs[10] = self.train_fish.GetFishPositionX(0) / 10
        # self.obs[11] = self.train_fish.GetFishPositionY(0) / 10
        # self.obs[12] = self.train_fish.GetFishVelocityX(0)
        # self.obs[13] = self.train_fish.GetFishVelocityY(0)
        # self.obs[12] = self.train_fish.GetFishPositionX(0)
        # self.obs[13] = self.train_fish.GetFishPositionY(0)
        # self.obs[14] = self.frequency / 10
        # self.obs[15] = self.lamb / 10
        self._get_obs = self.obs.astype(np.float32)
        print('train_obs: ', self._get_obs)

        return self._get_obs, {}

    def step(self, action):
        self.action_time_steps += 1
        self.frequency_new = action[0] + 3.0
        self.lamb_new = action[1] + 2.0
        print('train_fre: ', self.frequency_new)
        print('train_lam: ', self.lamb_new)

        file = open(f'action_{self.episode}.txt', 'a')
        file.write('action_time:  ')
        file.write(str(self.action_time))
        file.write('  ferquency:  ')
        file.write(str(self.frequency_new))
        file.write('  lamb:  ')
        file.write(str(self.lamb_new))
        file.write('\n')

        for i in range(20): # num_targets
            self.frequency = self.frequency + self.alpha * (self.frequency_new - self.frequency)
            self.lamb = self.lamb + self.alpha * (self.lamb_new - self.lamb)
            self.train_fish.SetFreq(self.frequency)
            self.train_fish.SetLambda(self.lamb)
            self.action_time += self.time_per_action / 20
            self.train_fish.RunCase(self.episode, self.action_time)

        self.fish_dx = self.train_fish.GetFishPositionX(0) - self.goal_x
        self.fish_dy = self.train_fish.GetFishPositionY(0) - self.goal_y
        print(' fish_dx: ', self.fish_dx, ' fish_dy: ', self.fish_dy)

        reward = 1 - np.sqrt(pow(self.fish_dx, 2) + pow(self.fish_dy, 2)) / np.sqrt(pow(self.fish_start_x - self.goal_x, 2) + pow(self.fish_start_y - self.goal_y, 2))
        print('reward: ', reward)
        reward_nan = checknan(reward)
        if reward_nan: # check reward nan
            reward = 0

        for i in range(32):
            self.obs[i] = self.train_fish.GetFishPressurePoint(i) / self.arr[i]  #
        # self.obs[10] = self.train_fish.GetFishPositionX(0) / 10
        # self.obs[11] = self.train_fish.GetFishPositionY(0) / 10
        # self.obs[12] = self.train_fish.GetFishVelocityX(0)
        # self.obs[13] = self.train_fish.GetFishVelocityY(0)
        # self.obs[12] = self.train_fish.GetFishPositionX(0)
        # self.obs[13] = self.train_fish.GetFishPositionY(0)
        # self.obs[14] = self.frequency / 10
        # self.obs[15] = self.lamb / 10
        self._get_obs = self.obs.astype(np.float32)

        done = False
        for i in range(32):
            if self.train_fish.GetFishPositionX(i) < 0 or self.train_fish.GetFishPositionX(i) > self.tank_lenth:
                done = True
                break
            if self.train_fish.GetFishPositionY(i) < 0 or self.train_fish.GetFishPositionY(i) > self.tank_height:
                done = True
                break
        if done == False:
            if self.action_time_steps > 99:
                done = True
            elif reward_nan:
                done = True
                reward = 0
            elif np.sqrt(pow(self.fish_dx, 2) + pow(self.fish_dy, 2)) <= 0.1:
                done = True
                reward = reward + 10.0
        # elif self.train_fish.GetFishPositionX(0) < 0.1 or self.train_fish.GetFishPositionY(0) > 0.55:
        #     done = True
        #     reward = reward - 10.0
        # elif self.train_fish.GetFishPositionY(0) < 0.05:
        #     done = True
        #     reward = reward - 5.0
        # else:
        #     done = False

        file = open(f'reward_{self.episode}.txt', 'a')
        file.write('action_time:  ')
        file.write(str(self.action_time))
        file.write('  reward:  ')
        file.write(str(reward))
        file.write('  fish_x  ')
        file.write(str(self.fish_dx + self.goal_x))
        file.write('  fish_y:  ')
        file.write(str(self.fish_dy + self.goal_y))
        file.write('\n')
        file.close()

        self.total_reward += reward

        if done == True:
            file = open('reward.txt', 'a')
            file.write('episode:  ')
            file.write(str(self.episode))
            file.write('  reward:  ')
            file.write(str(self.total_reward))
            file.write('\n')
            file.close()
            self.episode += 1

        return self._get_obs , reward, done, False, {}

    def render(self):
        return 0

    def _render_frame(self):
        return 0

    def close(self):
        return 0


if __name__=="__main__":
    env = FISHEnv()
    state_shape = env.observation_space.shape
    print('state_shape:', state_shape)
    action_shape = env.action_space.shape
    print('action_shape: ', action_shape)
    action_space = env.action_space
    print('action_space', action_space)
    init_state = env.reset()

    done = False
    while not done:
        act = env.action_space.sample()
        #print("this is sample action",act)
        s, r, done, _, _ = env.step(act)
        print("Action", act)
        print("Module State", s)
        print("Reward", r)
        print("Done", done)