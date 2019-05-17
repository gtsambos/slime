
import slime
import unittest
import re

class TestRecentHistoryGenerator(unittest.TestCase):

    scr = slime.RecentHistory()
    scr.add_event(1, 'early', 'start of simulation')
    scr.add_event(3, 'late', 'second event')
    scr.add_event(3, 'late', 'first event', start=True)
    scr.add_event(3, 'late', 'third event')
    scr.add_event(2, 'late', 'the middle of the script')
    scr.add_event(5, 'late', 'the end of the script')
    scr.add_event(4, 'early', 'more middle of the script')

    # """
    # Tests the methods in the RecentHistory class.
    # """

    def test_generation_order(self):
        scr = self.scr.dump_script()
        test_regex = re.compile("^\d+", flags=re.MULTILINE)
        generations = test_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        self.assertEqual(generations, sorted(generations))

    def test_time_order(self):
        scr = self.scr.dump_script()
        gen_regex = re.compile("^\d+", flags=re.MULTILINE)
        event_regex = re.compile(r"\d+ \bearly\b|\d+ \blate\b")
        generations = gen_regex.findall(scr)
        generations = [int(gen) for gen in generations]
        events = event_regex.findall(scr)
        for gen in generations:
            early_str = "%i early" % gen
            late_str = "%i late" % gen
            early_regex = re.compile(early_str)
            late_regex = re.compile(late_str)
            early = early_regex.findall(scr)
            late = late_regex.findall(scr)
            self.assertLess(len(early), 2)
            self.assertLess(len(late), 2)
            if len(early) == 1 and len(late) == 1:
                self.assertEqual(events.index(early) + 1, events.index(late))

    def test_ordering_inside_events(self):
        scr = self.scr.dump_script()
        event = '3 late(){'
        event_loc = scr.find(event)
        end_of_scr = scr[event_loc:]
        first_loc = end_of_scr.find('first event')
        second_loc = end_of_scr.find('second event')
        third_loc = end_of_scr.find('third event')
        self.assertLess(first_loc, second_loc)
        self.assertLess(second_loc, third_loc)

